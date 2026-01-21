#!/bin/csh -f
# 简洁版 GMT GRD 转 TIFF 脚本

# 参数检查
if ($#argv < 3) then
    echo "Usage: $0 input.grd color.cpt output.tif [scale_mode]"
    echo "scale_mode: pixel (1:1) or width like 15c, 10c (default: 15c)"
    exit 1
endif

set grd = $1
set cpt = $2
set out = $3

# 设置比例尺模式
if ($#argv == 4) then
    set mode = $4
else
    set mode = "15c"  # 默认15厘米宽度
endif

set base = $out:r

# 获取网格范围
set info = `gmt grdinfo $grd -C`
set xmin = $info[2]
set xmax = $info[3]
set ymin = $info[4]
set ymax = $info[5]
# echo $xmin, $xmax, $ymin, $ymax

# 设置投影
if ("$mode" == "pixel") then
    # 1:1像素模式 - 使用网格尺寸
    set nx = $info[10]
    set ny = $info[11]
    # set proj = "X${nx}/0"
    set proj = "X${nx}/Y${ny}"
    echo "Using 1:1 pixel mode (${nx}x${ny})"
else
    # 固定宽度模式
    set proj = "X${mode}"
    echo "Using fixed width: $mode"
endif

# 创建PS文件 (无边框)
echo "Converting to TIFF..."

if ("$mode" == "pixel") then
    gmt grdimage $grd -R$xmin/$xmax/$ymin/$ymax -J$proj -C$cpt  -Atest.tif=Gtiff
else
    #echo "gmt grdimage $grd -R$xmin/$xmax/$ymin/$ymax -J$proj -C$cpt # -B0"
    gmt grdimage $grd -R$xmin/$xmax/$ymin/$ymax -J$proj -C$cpt  -B0 > $base.ps
    gmt psconvert $base.ps -A -P -Tt
endif

# 重命名和清理
# mv $base.tiff $out
rm -f $base.ps gmt.history gmt.conf

echo "Done: $out created"
