时间步长Δt = 0.02， 网格间距为 x = 0.1(文件夹01)，0.05(文件夹02)，0.025(文件夹03)的均匀网格（即 Δx = Δy）;
网格间距Δx = Δy = 0.05， 时间步长Δt = 0.005(文件夹04)，0.01(文件夹05)，0.02(文件夹06)，0.04(文件夹07);
每个文件夹中：
dataout.csv为输出的解，纵坐标为j，横坐标为i，为一个时间步的解，然后按时间步堆叠在文件后面
errorout.csv为误差
exact.csv为解析解
adi.c为求解器
adi.py为曲线绘制后处理器
epsilon_1为时间1秒时的误差
eh.csv为随时间增长的误差

运行时cd到当前目录，然后generator回车（path中需要有gcc和python）

如果要单独运行求解器，则命令行输入adi [dx] [dt] [物理时间] [并行求解核心数量]

思考题：观察发现，误差随空间网格变小而变小，但是非线性，空间网格小到一定程度后，误差变小的速度放慢
误差随时间步长变小而增大，但是发现增大到一定程度后增幅变小直至趋于零
误差随物理时间先是线性减小，后来稳定到一定的值附近（可能与计算机双精度浮点数精度有关）
