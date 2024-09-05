# 导入需要的库
from pylab import meshgrid, cos
from numpy import c_, r_, diff, zeros, diag, matrix, identity, arange
from math import pi, sin
import netCDF4 as nc
import numpy as np
import pandas as pd
import hashlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import ListedColormap

d = nc.Dataset('E:\Data\GS1year\Temp\T_7.25.nc')
# lon = d0['lon'][:]
lat = d['lat'][:]
lon1 = d['lon'][:]
# d1中的经度是0.5，359.5]，而d2中的是[-80, 80]，因此将d1的经度做转换
lon2 = np.zeros_like(lon1)
for i in range(len(lon1)):
    if lon1[i] > 180:
        lon2[i] = lon1[i]-360  # 将西经[180.5, 359.5]转换为 -0.5 ~ -179.5
    else:
        lon2[i] = lon1[i]
lon0 = -lon2
lon = [round(num, 3) for num in lon0]
lon = np.array(lon)

l1 = 14  # z=-50m
l2 = 21  # z=-150m
l3 = 24  # z=-300m
l4 = 29  # z=-700m

t_i = d['t_i'][:]
t_H = d['t_truth'][:]
t_e = d['t_e'][:]

s_i = d['s_i'][:]
s_H = d['s_truth'][:]
s_e = d['s_e'][:]

colormap = 'coolwarm'   # 'RdYlBu' 'Spectral' 'coolwarm'

''' 画图——温度对比图（isQG、HYCOM、eSQG、truth） 0、100、200、300m'''

fig1, ax1 = plt.subplots(3, 4, figsize=(13, 9))

levels1 = np.arange(14, 26, 0.5)
levels12 = np.arange(35.2, 37.2, 0.1)
ax1[0, 0].contourf(abs(lon), abs(lat), t_i[l1, :, :], levels1, cmap=colormap, extend='both')
contourf1 = ax1[0, 1].contourf(abs(lon), abs(lat), t_H[l1, :, :], levels1, cmap=colormap, extend='both')
ax1[0, 2].contourf(abs(lon), abs(lat), s_i[l1, :, :], levels12, cmap=colormap, extend='both')  # 'autumn'
contour1 = ax1[0, 2].contour(abs(lon), abs(lat), t_H[l1, :, :], colors='w', linewidths=0.7)
contourf12 = ax1[0, 3].contourf(abs(lon), abs(lat), s_H[l1, :, :], levels12, cmap=colormap, extend='both')

ax1[0, 0].invert_xaxis()  # 反转坐标轴
ax1[0, 1].invert_xaxis()  # 反转坐标轴
ax1[0, 2].invert_xaxis()  # 反转坐标轴
ax1[0, 3].invert_xaxis()  # 反转坐标轴
ax1[0, 0].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [0, 0]窗口
ax1[0, 0].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))
ax1[0, 1].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))
ax1[0, 1].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))
ax1[0, 2].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))
ax1[0, 2].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))
ax1[0, 3].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))
ax1[0, 3].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))

levels2 = np.arange(11, 21, 0.5)
levels22 = np.arange(35.2, 37.1, 0.1)
ax1[1, 0].contourf(abs(lon), abs(lat), t_i[l2, :, :], levels2, cmap=colormap, extend='both')
contourf2 = ax1[1, 1].contourf(abs(lon), abs(lat), t_H[l2, :, :], levels2, cmap=colormap, extend='both')
ax1[1, 2].contourf(abs(lon), abs(lat), s_i[l2, :, :], levels22, cmap=colormap, extend='both')  # 'autumn'
contour2 = ax1[1, 2].contour(abs(lon), abs(lat), t_H[l2, :, :], colors='w', linewidths=0.7)
contourf22 = ax1[1, 3].contourf(abs(lon), abs(lat), s_H[l2, :, :], levels22, cmap=colormap, extend='both')

ax1[1, 0].invert_xaxis()  # 反转坐标轴
ax1[1, 1].invert_xaxis()  # 反转坐标轴
ax1[1, 2].invert_xaxis()  # 反转坐标轴
ax1[1, 3].invert_xaxis()  # 反转坐标轴

ax1[1, 0].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [0, 0]窗口
ax1[1, 0].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))
ax1[1, 1].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))
ax1[1, 1].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))
ax1[1, 2].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))
ax1[1, 2].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))
ax1[1, 3].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))
ax1[1, 3].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))

levels3 = np.arange(10, 20, 0.5)
levels32 = np.arange(35, 37., 0.1)
ax1[2, 0].contourf(abs(lon), abs(lat), t_i[l3, :, :], levels3, cmap=colormap, extend='both')
contourf3 = ax1[2, 1].contourf(abs(lon), abs(lat), t_H[l3, :, :], levels3, cmap=colormap, extend='both')
ax1[2, 2].contourf(abs(lon), abs(lat), s_i[l3, :, :], levels32, cmap=colormap, extend='both')  # 'autumn'
contour3 = ax1[2, 2].contour(abs(lon), abs(lat), t_H[l3, :, :], colors='w', linewidths=0.7)
contourf32 = ax1[2, 3].contourf(abs(lon), abs(lat), s_H[l3, :, :], levels32, cmap=colormap, extend='both')

ax1[2, 0].invert_xaxis()  # 反转坐标轴
ax1[2, 1].invert_xaxis()  # 反转坐标轴
ax1[2, 2].invert_xaxis()  # 反转坐标轴
ax1[2, 3].invert_xaxis()  # 反转坐标轴

ax1[2, 0].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [0, 0]窗口
ax1[2, 0].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))
ax1[2, 1].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))
ax1[2, 1].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))
ax1[2, 2].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))
ax1[2, 2].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))
ax1[2, 3].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))
ax1[2, 3].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))

cax1 = fig1.add_axes([0.15, 0.63, 0.3, 0.015])  # 创建颜色条的位置和大小
bar1 = plt.colorbar(contourf1, cax=cax1, orientation='horizontal')
bar1.set_ticks([16, 18, 20, 22, 24])  # 设置颜色条刻度

cax12 = fig1.add_axes([0.57, 0.63, 0.3, 0.015])  # 创建颜色条的位置和大小
bar12 = plt.colorbar(contourf12, cax=cax12, orientation='horizontal')
bar12.set_ticks([35.4, 35.7, 36, 36.3, 36.6, 36.9])  # 设置颜色条刻度

cax2 = fig1.add_axes([0.15, 0.34, 0.3, 0.015])  # 创建颜色条的位置和大小
bar2 = plt.colorbar(contourf2, cax=cax2, orientation='horizontal')
bar2.set_ticks([12, 14, 16, 18, 20])  # 设置颜色条刻度

cax22 = fig1.add_axes([0.57, 0.34, 0.3, 0.015])  # 创建颜色条的位置和大小
bar22 = plt.colorbar(contourf22, cax=cax22, orientation='horizontal')
bar22.set_ticks([35.3, 35.6, 35.9, 36.2, 36.5, 36.8])  # 设置颜色条刻度

cax3 = fig1.add_axes([0.15, 0.05, 0.3, 0.015])  # 创建颜色条的位置和大小
bar3 = plt.colorbar(contourf3, cax=cax3, orientation='horizontal')
bar3.set_ticks([11, 13, 15, 17, 19])  # 设置颜色条刻度

cax32 = fig1.add_axes([0.57, 0.05, 0.3, 0.015])  # 创建颜色条的位置和大小
bar32 = plt.colorbar(contourf32, cax=cax32, orientation='horizontal')
bar32.set_ticks([35.2, 35.5, 35.8, 36.1, 36.4, 36.7])  # 设置颜色条刻度

ax1[0, 0].set_yticks([31, 34, 37, 40])  # 设置y轴刻度标签
ax1[0, 1].set_yticks([31, 34, 37, 40])
ax1[0, 2].set_yticks([31, 34, 37, 40])
ax1[0, 3].set_yticks([31, 34, 37, 40])
ax1[1, 0].set_yticks([31, 34, 37, 40])
ax1[1, 1].set_yticks([31, 34, 37, 40])
ax1[1, 2].set_yticks([31, 34, 37, 40])
ax1[1, 3].set_yticks([31, 34, 37, 40])
ax1[2, 1].set_yticks([31, 34, 37, 40])
ax1[2, 0].set_yticks([31, 34, 37, 40])
ax1[2, 2].set_yticks([31, 34, 37, 40])
ax1[2, 3].set_yticks([31, 34, 37, 40])

plt.subplots_adjust(wspace=0.35, hspace=0.55)

# 添加整个窗口的 x 轴和 y 轴标签
fig1.text(0.5, 0.03, 'Longitude', ha='center', font={'family': 'Arial', 'size': 16})
fig1.text(0.02, 0.5, 'Latitude', va='center', rotation='vertical', font={'family': 'Arial', 'size': 16})

fig1.text(0.2, 0.90, '$Temp^{isQG}$', ha='center', font={'family': 'Arial', 'size': 14})
fig1.text(0.42, 0.90, '$Temp^{HYCOM}$', ha='center', font={'family': 'Arial', 'size': 14})
fig1.text(0.62, 0.90, '$Salinity^{isQG}$', ha='center', font={'family': 'Arial', 'size': 14})
fig1.text(0.83, 0.90, '$Salinity^{HYCOM}$', ha='center', font={'family': 'Arial', 'size': 14})

fig1.text(0.06, 0.76, '150m', va='center', rotation='vertical', font={'family': 'Arial', 'size': 14})
fig1.text(0.06, 0.5, '300m', va='center', rotation='vertical', font={'family': 'Arial', 'size': 14})
fig1.text(0.06, 0.2, '700m', va='center', rotation='vertical', font={'family': 'Arial', 'size': 14})

# ********************** profile
z = d['depth'][:]
x = 98  # 52.16°W
y = 88  # 33.52°N

ti1 = np.mean(t_i, axis=(1, 2))
te1 = np.mean(t_e, axis=(1, 2))
tH1 = np.mean(t_H, axis=(1, 2))

delt_ti = diff(t_i, axis=0)
delt_te = diff(t_e, axis=0)
delt_tH = diff(t_H, axis=0)

dz = diff(-z)

gradient_i = delt_ti/dz.reshape(-1, 1, 1)
gradient_e = delt_te/dz.reshape(-1, 1, 1)
gradient_H = delt_tH/dz.reshape(-1, 1, 1)


''' 画图——温度对比图（isQG、HYCOM、eSQG、truth） 0、100、200、300m'''

fig1, ax1 = plt.subplots(2, 2, figsize=(9, 9))

levels1 = np.arange(4, 25, 1.4)
levels2 = np.arange(33.5, 39.9, 0.4)

contourf1 = ax1[0, 0].contourf(abs(lon), -z[0: 35], t_i[:35, y, :], levels1, cmap=colormap1, extend='both')
ax1[0, 0].set_title('(a)  Temperature (isQG)', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[0, 0].invert_xaxis()  # 反转坐标轴
ax1[0, 0].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f m'))  # 给纬度加单位  [0, 0]窗口
ax1[0, 0].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))

contourf2 = ax1[0, 1].contourf(abs(lon), -z[0: 35], s_i[:35, y, :], levels2, cmap=colormap, extend='both')
ax1[0, 1].set_title('(b)  Salinity (isQG)', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[0, 1].invert_xaxis()  # 反转坐标轴
ax1[0, 1].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f m'))  # 给纬度加单位  [0, 0]窗口
ax1[0, 1].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))


contourf11 = ax1[1, 0].contourf(abs(lon), -z[0: 35], t_H[:35, y, :], levels1, cmap=colormap1, extend='both')
ax1[1, 0].set_title('(c)  Temperature (HYCOM)', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[1, 0].invert_xaxis()  # 反转坐标轴
ax1[1, 0].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f m'))  # 给纬度加单位  [0, 0]窗口
ax1[1, 0].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))


contourf21 = ax1[1, 1].contourf(abs(lon), -z[0: 35], s_H[:35, y, :], levels2, cmap=colormap, extend='both')
ax1[1, 1].set_title('(d)  Salinity (HYCOM)', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[1, 1].invert_xaxis()  # 反转坐标轴
ax1[1, 1].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f m'))  # 给纬度加单位  [0, 0]窗口
ax1[1, 1].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))


cax1 = fig1.add_axes([0.13, 0.05, 0.32, 0.015])  # 创建颜色条的位置和大小
bar1 = plt.colorbar(contourf11, cax=cax1, orientation='horizontal')  # vertical
bar1.set_ticks([7, 10, 13, 16, 19, 22])  # 设置颜色条刻度

cax2 = fig1.add_axes([0.57, 0.05, 0.32, 0.015])  # 创建颜色条的位置和大小[左，底，宽，高]
bar2 = plt.colorbar(contourf21, cax=cax2, orientation='horizontal')  # vertical
bar2.set_ticks([34, 35, 36, 37, 38, 39])  # 设置颜色条刻度
bar2.ax.tick_params(labelsize=10)

plt.subplots_adjust(wspace=0.3)

# 添加整个窗口的 x 轴和 y 轴标签
fig1.text(0.5, 0.01, 'Longitude', ha='center', font={'family': 'Arial', 'size': 16})
fig1.text(0.02, 0.5, 'Depth', va='center', rotation='vertical', font={'family': 'Arial', 'size': 16})
fig1.text(0.5, 0.94, 'July', ha='center', font={'family': 'Arial', 'size': 18})

# *********** correlation

At = np.corrcoef(t_i[:, y, :], t_H[:, y, :])
As = np.corrcoef(s_i[:, y, :], s_H[:, y, :])
cor_t = np.zeros(len(z))
cor_s = np.zeros(len(z))
for k in range(len(z)):
    cor_t[k] = At[k, 38+k]
    cor_s[k] = As[k, 38+k]

def R2(y_obs, y_sim):
    y_obs_mean = np.mean(y_obs)
    y_sim_mean = np.mean(y_sim)
    sum1 = 0
    sum2 = 0
    sum3 = 0
    for i in range(len(y_obs)):
        sum1 = sum1 + (y_sim[i] - y_sim_mean) * (y_obs[i] - y_obs_mean)
        sum2 = sum2 + ((y_sim[i] - y_sim_mean) ** 2)
        sum3 = sum3 + ((y_obs[i] - y_obs_mean) ** 2)
    R2 = (sum1 / ((sum2 ** 0.5) * (sum3 ** 0.5))) ** 2
    return R2


def get_mse(records_real, records_predict):
    """
    均方误差 估计值与真值 偏差
    """
    return sum([(x - y) ** 2 for x, y in zip(records_real, records_predict)]) / len(records_real)


def get_rmse(records_real, records_predict):
    """
    均方根误差：是均方误差的算术平方根
    """
    mse = get_mse(records_real, records_predict)
    return np.sqrt(mse)


R2_t = np.zeros(len(z))
R2_s = np.zeros(len(z))
RMSE_t = np.zeros(len(z))
RMSE_s = np.zeros(len(z))

for k in range(len(z)):
    R2_t[k] = R2(t_H[k, y, :], t_i[k, y, :])
    R2_s[k] = R2(s_H[k, y, :], s_i[k, y, :])
    RMSE_t[k] = get_rmse(t_H[k, y, :], t_i[k, y, :])
    RMSE_s[k] = get_rmse(s_H[k, y, :], s_i[k, y, :])

fig1, ax = plt.subplots(1, 2, figsize=(8, 4))
ax[0].plot(cor_t, -z, color='orchid', label='Temperature')  # 画垂直分层, lightcoral
ax[0].plot(cor_s, -z, color='cornflowerblue', label='Salinity')  # royalblue
ax[0].legend(loc='best')  # 添加图例
# ax[0].set_xlim(0.6, 1)  # Y轴范围2-7
ax[0].set_title(" Correlation of temperature", font={'family': 'Arial', 'size': 16})  # 添加轴标签

ax[1].plot(RMSE_t, -z, color='lightsalmon', label='Temperature')  # 画垂直分层, lightcoral
ax[1].plot(RMSE_s, -z, color='slateblue', label='Salinity')  # rosybrown
ax[1].legend(loc='best')  # 添加图例
# x[1].set_xlim(0, 0)  # Y轴范围2-7
ax[1].set_title(" Root-mean-square error", font={'family': 'Arial', 'size': 16})  # 添加轴标签
fig1.text(0.02, 0.5, 'Depth(m)', va='center', rotation='vertical', font={'family': 'Arial', 'size': 16})

plt.show()
