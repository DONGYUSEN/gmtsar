#!/bin/bash
# GMT绘制笛卡尔矢量（起点+终点）+ 相关性气泡图（修复-R参数错误版）
# 输入：freq_xcorr.dat（列1=X, 列2=deltaX, 列3=Y, 列4=deltaY, 列5=相关性）
# 输出：freq_xcorr.pdf（自动生成）

# ===================== 1. 基础配置 =====================
INPUT_DATA="freq_xcorr.dat"
OUTPUT_PDF="${INPUT_DATA%.dat}.pdf"
TMP_PS="temp_plot.ps"
gmt set PS_MEDIA=A0
# 检查输入文件
if [ ! -f "$INPUT_DATA" ]; then
    echo "错误：输入文件 $INPUT_DATA 不存在！"
    exit 1
fi

# ===================== 2. 数据预处理（核心：修复-R格式） =====================
# 步骤1：提取X/Y的最大值、最小值
X_MAX=$(awk 'BEGIN{max_x=0} {if($1!="" && $1+0>max_x) max_x=$1+0} END{print max_x}' $INPUT_DATA)
Y_MAX=$(awk 'BEGIN{max_y=0} {if($3!="" && $3+0>max_y) max_y=$3+0} END{print max_y}' $INPUT_DATA)
X_MIN=$(awk 'BEGIN{min_x=1000000} {if($1!="" && $1+0<min_x) min_x=$1+0} END{print min_x}' $INPUT_DATA)
Y_MIN=$(awk 'BEGIN{min_y=1000000} {if($3!="" && $3+0<min_y) min_y=$3+0} END{print min_y}' $INPUT_DATA)

# 步骤2：计算REGION范围（严格遵守west/east/south/north，且west<east、south<north）
X_REGION_MIN=$((X_MIN - 500))  # X起始：X_MIN-500
X_REGION_MAX=$((X_MAX + 500))  # X结束：X_MAX+500
Y_REGION_MIN=$((Y_MIN - 500))  # Y起始（south）：Y_MIN-500
Y_REGION_MAX=$((Y_MAX + 1000))  # Y结束（north）：Y_MAX+500
# 关键修复：REGION格式严格为 west/east/south/north（X_MIN-500/X_MAX+500/Y_MIN-500/Y_MAX+500）
REGION="$X_REGION_MIN/${X_REGION_MAX}/${Y_REGION_MIN}/${Y_REGION_MAX}"

# 步骤4：计算图幅纵横比（匹配Y/X范围）
BASE_WIDTH=16  # 基准宽度（cm）
X_RANGE=$((X_REGION_MAX - X_REGION_MIN + 500))
Y_RANGE=$((Y_REGION_MAX - Y_REGION_MIN + 1000))
ASPECT_RATIO=$(echo "scale=6; $Y_RANGE / $X_RANGE" | bc)
PROJ_HEIGHT=$(echo "$BASE_WIDTH * $ASPECT_RATIO" | bc)
PROJ="X+${BASE_WIDTH}c/-${PROJ_HEIGHT}c"  # 补全-J前缀，GMT必选

# 步骤5：偏移量缩放系数（保留原有逻辑）
DELTAX_MAX=$(awk 'BEGIN{max=0} {abs=($2>=0)?$2:-$2; if(abs>max) max=abs} END{print max}' $INPUT_DATA)
DELTAY_MAX=$(awk 'BEGIN{max=0} {abs=($4>=0)?$4:-$4; if(abs>max) max=abs} END{print max}' $INPUT_DATA)
DELTA_MAX=$(echo "$DELTAX_MAX $DELTAY_MAX" | awk '{if($1>$2) print $1; else print $2}')
TARGET_ARROW_LENGTH=0.6
SCALE_FACTOR=1

# 步骤6：转换为笛卡尔矢量格式
awk -v scale=$SCALE_FACTOR '
    {
        start_x = $1;
        start_y = $3;
        end_x   = start_x - $2;
        end_y   = start_y - $4;
        print start_x, start_y, end_x, end_y, $5;
    }
' $INPUT_DATA > vector_data.dat
awk '{print $1,$2,$3,$4}' vector_data.dat > vector_data2.dat

# 步骤7：提取相关性范围
CORR_MIN=$(awk 'BEGIN{min_c=1000} {if($5!="" && $5+0<min_c) min_c=$5+0} END{print min_c}' vector_data.dat)
CORR_MAX=$(awk 'BEGIN{max_c=0} {if($5!="" && $5+0>max_c) max_c=$5+0} END{print max_c}' vector_data.dat)
CORR_MIN=${CORR_MIN:-40}
CORR_MAX=${CORR_MAX:-70}

# ===================== 3. 绘图参数 =====================
ARROW_STYLE="-SV0.15c+ea+s -W0.5p -Gred "

# 生成相关性气泡大小文件
awk -v min=$CORR_MIN -v max=$CORR_MAX '
    BEGIN{scale=(0.55-0.15)/(max-min)}
    {print $1,$2,0.15+($5-min)*scale}
' vector_data.dat > corr_size.tmp

# ===================== 4. GMT绘图（修复-R错误） =====================
# 初始化绘图
gmt psxy -R$REGION -J$PROJ -T -K > $TMP_PS

# 绘制底图（仅保留图框+外部标注，无网格线）
# 关键：移除fg（网格线），仅保留刻度标注，-BWSne保留外框
gmt psbasemap -R$REGION -J$PROJ -Bxa1000+l"X" -Bya2000+l"Y" -BWSne -K -O >> $TMP_PS

# 绘制笛卡尔矢量
gmt psxy vector_data2.dat -R$REGION -J$PROJ $ARROW_STYLE -K -O >> $TMP_PS

# 绘制相关性气泡
gmt psxy corr_size.tmp -R$REGION -J$PROJ -Sc -W0.15p,black -Gblue@70 -K -O >> $TMP_PS

# 添加标题（适配反转后的Y轴）
TITLE_X=$(echo "scale=2; ($X_REGION_MIN + $X_REGION_MAX)/scale" | bc)
TITLE_Y=$(echo "$Y_REGION_MAX + ($Y_RANGE * 0.05)" | bc)  # 标题在Y轴顶部
gmt pstext -R$REGION -J$PROJ -N -K -O -F+f14p,Helvetica,black+jTC >> $TMP_PS << EOF
$TITLE_X $TITLE_Y offset+corr
EOF

# 添加图例
#LEGEND_X=$(echo "$X_REGION_MIN + ($X_RANGE * 0.05)" | bc)
#LEGEND_Y=$(echo "$Y_REGION_MAX - ($Y_RANGE * 0.1)" | bc)
#gmt pstext -R$REGION -J$PROJ -N -K -O -F+f12p,Helvetica,black+jBL >> $TMP_PS << EOF
#$LEGEND_X $LEGEND_Y Corr：${CORR_MIN}~ ${CORR_MAX}
#EOF

# 结束绘图
gmt psxy -R$REGION -J$PROJ -T -O >> $TMP_PS

# ===================== 5. 转换PDF+清理 =====================
gmt psconvert $TMP_PS -Tf -A -P -E300 -F$OUTPUT_PDF

# 恢复GMT默认设置（避免影响后续绘图）

# 清理临时文件
rm -f $TMP_PS corr_size.tmp vector_data.dat vector_data2.dat gmt.history gmt.conf gmt.command
gmt set PS_MEDIA=A4
# ===================== 6. 验证输出 =====================
echo "===================== 运行完成（修复-R错误） ====================="
echo "✅ 核心修复："
echo "   1. REGION格式：$REGION（严格遵守west/east/south/north，west<east、south<north）"
echo "   2. 图幅纵横比：${BASE_WIDTH}cm × ${PROJ_HEIGHT}cm（匹配数据范围）"
echo "✅ 输出文件：$OUTPUT_PDF"
echo "======================================================================"
