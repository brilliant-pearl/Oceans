# 导入需要的库

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

d15 = nc.Dataset('path\T_1.5.nc')
d115 = nc.Dataset('path\T_1.15.nc')
d125 = nc.Dataset('path\T_1.25.nc')
d25 = nc.Dataset('path\T_2.5.nc')
d215 = nc.Dataset('path\T_2.15.nc')
d225 = nc.Dataset('path\T_2.25.nc')
d35 = nc.Dataset('path\T_3.5.nc')
d315 = nc.Dataset('path\T_3.15.nc')
d325 = nc.Dataset('path\T_3.25.nc')
d45 = nc.Dataset('path\T_4.5.nc')
d415 = nc.Dataset('path\T_4.15.nc')
d425 = nc.Dataset('path\T_4.25.nc')
d55 = nc.Dataset('path\T_5.5.nc')
d515 = nc.Dataset('path\T_5.15.nc')
d525 = nc.Dataset('path\T_5.25.nc')
d65 = nc.Dataset('path\T_6.5.nc')
d615 = nc.Dataset('path\T_6.15.nc')
d625 = nc.Dataset('path\T_6.25.nc')
d75 = nc.Dataset('path\T_7.5.nc')
d715 = nc.Dataset('path\T_7.15.nc')
d725 = nc.Dataset('path\T_7.25.nc')
d85 = nc.Dataset('path\T_8.5.nc')
d815 = nc.Dataset('path\T_8.15.nc')
d825 = nc.Dataset('path\T_8.25.nc')
d95 = nc.Dataset('path\T_9.5.nc')
d915 = nc.Dataset('path\T_9.15.nc')
d925 = nc.Dataset('path\T_9.25.nc')
d105 = nc.Dataset('path\T_10.5.nc')
d1015 = nc.Dataset('path\T_10.15.nc')
d1025 = nc.Dataset('path\T_10.25.nc')

d1105 = nc.Dataset('path\T_11.5.nc')
d1115 = nc.Dataset('path\T_11.15.nc')
d1125 = nc.Dataset('path\T_11.25.nc')

d1205 = nc.Dataset('path\T_12.5.nc')
d1215 = nc.Dataset('path\T_12.15.nc')
d1225 = nc.Dataset('path\T_12.25.nc')

# lon = d0['lon'][:]
z = d15['depth'][:]
lat = d15['lat'][:]
lon1 = d15['lon'][:]

lon2 = np.zeros_like(lon1)
for i in range(len(lon1)):
    if lon1[i] > 180:
        lon2[i] = lon1[i]-360  # 将西经[180.5, 359.5]转换为 -0.5 ~ -179.5
    else:
        lon2[i] = lon1[i]
lon0 = -lon2
lon = [round(num, 3) for num in lon0]
lon = np.array(lon)

T = np.empty((4, len(z), len(lat), len(lon)))
S = np.empty((4, len(z), len(lat), len(lon)))


T[0, :, :, :] = d115['t_i'][:]
S[0, :, :, :] = d115['s_i'][:]

T[1, :, :, :] = d415['t_i'][:]
S[1, :, :, :] = d415['s_i'][:]

T[2, :, :, :] = d715['t_i'][:]
S[2, :, :, :] = d715['s_i'][:]

T[3, :, :, :] = d1015['t_i'][:]
S[3, :, :, :] = d1015['s_i'][:]


TH = np.empty((4, len(z), len(lat), len(lon)))
SH = np.empty((4, len(z), len(lat), len(lon)))

TH[0, :, :, :] = d115['t_truth'][:]
SH[0, :, :, :] = d115['s_truth'][:]

TH[1, :, :, :] = d415['t_truth'][:]
SH[1, :, :, :] = d415['s_truth'][:]

TH[2, :, :, :] = d715['t_truth'][:]
SH[2, :, :, :] = d715['s_truth'][:]

TH[3, :, :, :] = d1015['t_truth'][:]
SH[3, :, :, :] = d1015['s_truth'][:]


x = 98  # 52.16°W
y = 88  # 33.52°N

t = np.mean(T, axis=(2, 3))
s = np.mean(S, axis=(2, 3))
th = np.mean(TH, axis=(2, 3))
sh = np.mean(SH, axis=(2, 3))



''' 画图——温度对比图（isQG、HYCOM、eSQG、truth） 0、100、200、300m'''

fig1, ax = plt.subplots(2, 4, figsize=(11, 8))

ax[0, 0].plot(t[0, 0:35], -z[0:35], color='plum', marker='.', label='isQG')
ax[0, 0].plot(th[0, 0:35], -z[0:35], color='slategrey', label='HYCOM')
ax[0, 0].set_title('(a) Temperature (Jan)', font={'family': 'Arial', 'size': 14})
ax[0, 0].legend(loc='best')
ax[0, 1].plot(t[1, 0:35], -z[0:35], color='plum', marker='.', label='isQG')
ax[0, 1].plot(th[1, 0:35], -z[0:35], color='slategrey', label='HYCOM')
ax[0, 1].set_title('(b) Temperature (Apr)', font={'family': 'Arial', 'size': 14})
ax[0, 1].legend(loc='best')
ax[0, 2].plot(t[2, 0:35], -z[0:35], color='plum', marker='.', label='isQG')
ax[0, 2].plot(th[2, 0:35], -z[0:35], color='slategrey', label='HYCOM')
ax[0, 2].set_title('(c) Temperature (Jul)', font={'family': 'Arial', 'size': 14})
ax[0, 2].legend(loc='best')
ax[0, 3].plot(t[2, 0:35], -z[0:35], color='plum', marker='.', label='isQG')
ax[0, 3].plot(th[2, 0:35], -z[0:35], color='slategrey', label='HYCOM')
ax[0, 3].set_title('(d) Temperature (Oct)', font={'family': 'Arial', 'size': 14})
ax[0, 3].legend(loc='best')
ax[1, 0].plot(s[0, 0:35], -z[0:35], color='cornflowerblue', marker='.', label='isQG')
ax[1, 0].plot(sh[0, 0:35], -z[0:35], color='slategrey', label='HYCOM')
ax[1, 0].set_title('(e) Salinity (Jan)', font={'family': 'Arial', 'size': 14})
ax[1, 0].legend(loc='best')
ax[1, 1].plot(s[1, 0:35], -z[0:35], color='cornflowerblue', marker='.', label='isQG')
ax[1, 1].plot(sh[1, 0:35], -z[0:35], color='slategrey', label='HYCOM')
ax[1, 1].set_title('(f) Salinity (Apr)', font={'family': 'Arial', 'size': 14})
ax[1, 1].legend(loc='best')
ax[1, 2].plot(s[2, 0:35], -z[0:35], color='cornflowerblue', marker='.', label='isQG')
ax[1, 2].plot(sh[2, 0:35], -z[0:35], color='slategrey', label='HYCOM')
ax[1, 2].set_title('(g) Salinity (Jul)', font={'family': 'Arial', 'size': 14})
ax[1, 2].legend(loc='best')
ax[1, 3].plot(s[2, 0:35], -z[0:35], color='cornflowerblue', marker='.', label='isQG')
ax[1, 3].plot(sh[2, 0:35], -z[0:35], color='slategrey', label='HYCOM')
ax[1, 3].set_title('(h) Salinity (Oct)', font={'family': 'Arial', 'size': 14})
ax[1, 3].legend(loc='best')
fig1.text(0.04, 0.5, 'Depth(m)', va='center', rotation='vertical', font={'family': 'Arial', 'size': 14})

plt.subplots_adjust(wspace=0.35, hspace=0.3)

l1 = 14  # z=-50m
l2 = 21  # z=-150m
l3 = 24  # z=-300m
l4 = 29  # z=-700m

t = np.mean(T, axis=(2, 3))
s = np.mean(S, axis=(2, 3))

u = np.arange(1, 37)

# ********** one year

fig1, ax = plt.subplots(2, 1, figsize=(7, 8))

ax[0].plot(u, t[:, 0], color='mediumpurple', marker='+', label='z= 0m')
ax[0].plot(u, t[:, l1], color='slategrey', marker='.', label='z=-50m')
ax[0].plot(u, t[:, l2], color='darkseagreen', marker='1', label='z=-150m')
ax[0].plot(u, t[:, l3], color='cornflowerblue', marker='+', label='z=-300m')
ax[0].plot(u, t[:, l4], color='plum', marker='*', label='z=-700m')  # burlywood
ax[0].set_title('(a) Temperature', font={'family': 'Arial', 'size': 16})
ax[0].set_xticks([2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35], ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
ax[0].set_yticks([11, 14, 17, 20, 23, 26])
ax[0].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°C'))
# ax[0].set_ylabel('(deg  C)', labelpad=-36, y=1.02, rotation=0, font={'family': 'Arial', 'size': 9})  # 添加标题
ax[0].tick_params(axis="both", which="major", direction="in", width=1, length=5, pad=5)
ax[0].legend(loc='best', prop={'size': 8})
ax[0].grid(color='lightgray', linewidth=0.5, alpha=0.8)

ax[1].plot(u, s[:, 0], color='mediumpurple', marker='+', label='z=0m')
ax[1].plot(u, s[:, l1], color='slategrey', marker='.', label='z=-50m')
ax[1].plot(u, s[:, l2], color='darkseagreen', marker='1', label='z=-150m')
ax[1].plot(u, s[:, l3], color='cornflowerblue', marker='+', label='z=-300m')
ax[1].plot(u, s[:, l4], color='plum', marker='*', label='z=-700m')
ax[1].set_title('(b) Salinity', font={'family': 'Arial', 'size': 16})
ax[1].set_xticks([2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35], ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
ax[1].legend(loc='best')
ax[1].grid(color='lightgray', linewidth=0.5, alpha=0.8)

plt.show()
