import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import copy
import matplotlib.ticker as mticker

d1 = nc.Dataset('path_isQG.nc')
d2 = nc.Dataset('path_eSQG.nc')
d3 = nc.Dataset('path_GS.nc')
lat = d1['lat'][:]
lon = d3['lon'][:]
z = d1['depth'][:]

u_i = d1['u'][:]
u_e = d2['u'][:]
u_H = d3['u'][:]

v_i = d1['v'][:]
v_e = d2['v'][:]
v_H = d3['v'][:]

rho_i = d1['rhot'][:]
rho_e = d2['density'][:]
rho_H = d3['rho'][:]

rho_ii = d1['rhoi'][:]
rho_is = d1['rhos'][:]
l = 25  # 画图要展示的深度

levels1 = np.arange(-0.4, 0.48, 0.05)
levels2 = np.arange(-0.2, 0.2, 0.01)
fig1, ax1 = plt.subplots(2, 2, figsize=(8.2, 7))
contourf11 = ax1[0, 0].contourf(abs(lon), abs(lat), rho_is[l, :, :], levels1, cmap='coolwarm', extend='both')
contourf12 = ax1[0, 1].contourf(abs(lon), abs(lat), rho_ii[l, :, :], levels1, cmap='coolwarm', extend='both')
contourf13 = ax1[1, 0].contourf(abs(lon), abs(lat), rho_i[l, :, :], levels1, cmap='coolwarm', extend='both')
contourf14 = ax1[1, 1].contourf(abs(lon), abs(lat), rho_H[l, :, :], levels1, cmap='coolwarm', extend='both')
ax1[0, 0].set_title('(a)  $ \\rho^{s}_{z=-150}$', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[0, 1].set_title('(b)  $ \\rho^{i}_{z=-150}$', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[1, 0].set_title('(c)  $ \\rho^{isQG}_{z=-150}$', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[1, 1].set_title('(d)  $ \\rho^{HYCOM}_{z=-150}$', font={'family': 'Arial', 'size': 12})  # 添加标题

x = np.zeros_like(lon)
y = x+33.52
ax1[1, 0].plot(abs(lon), y, color='w', linestyle='--')  # 在lat=33.52°N位置处加一条白色虚线作为标记
ax1[1, 1].plot(abs(lon), y, color='w', linestyle='--')

ax1[0, 0].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [0, 0]窗口
ax1[0, 0].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位
ax1[0, 1].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [0, 1]
ax1[0, 1].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位
ax1[1, 0].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [1, 0]
ax1[1, 0].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位
ax1[1, 1].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [1, 1]
ax1[1, 1].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位

ax1[0, 0].invert_xaxis()  # 反转坐标轴
ax1[0, 1].invert_xaxis()  # 反转坐标轴
ax1[1, 0].invert_xaxis()  # 反转坐标轴
ax1[1, 1].invert_xaxis()  # 反转坐标轴

ax1[0, 0].set_yticks([31, 34, 37, 40])  # 设置y轴刻度标签
ax1[0, 1].set_yticks([31, 34, 37, 40])
ax1[1, 0].set_yticks([31, 34, 37, 40])
ax1[1, 1].set_yticks([31, 34, 37, 40])

cax = fig1.add_axes([0.92, 0.2, 0.02, 0.6])  # 创建颜色条的位置和大小
bar = plt.colorbar(contourf14, cax=cax, orientation='vertical')
bar.set_ticks([-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3])  # 设置颜色条刻度

fig1.text(0.5, 0.02, 'Longitude', ha='center', font={'family': 'Arial', 'size': 16})
fig1.text(0.02, 0.5, 'Latitude', va='center', rotation='vertical', font={'family': 'Arial', 'size': 16})

# *********** rho+uv

l1 = 24
l2 = 29  # 画图要展示的深度

lon11 = d1['lon'][:]

lon2 = np.zeros_like(lon11)
for i in range(len(lon11)):
    if lon11[i] > 180:
        lon2[i] = lon11[i]-360
    else:
        lon2[i] = lon11[i]
lon = -lon2

lon2 = np.float32(np.arange(44, 60.1, 0.4))
lat2 = np.float32(np.arange(30, 42.1, 0.4))
lon3 = [round(num, 3) for num in lon2]
lat1 = [round(num, 3) for num in lat2]
lon1 = list(reversed(lon3))
lat0 = copy.deepcopy(lat)
lon0 = copy.deepcopy(lon)
lon0 = [round(num, 3) for num in lon0]
lat0 = [round(num, 3) for num in lat0]


def filter3D_lonlat(d, lon0, lat0, lon, lat):
    import numpy as np
    '''
    :param data: A array of 2D or 3D
    :param lon0: A 1D array of longitude value which is corresponded to d.
    :param lat0: A 1D array of latitude value which is corresponded to d.
    :param lon: A 1D array of longitude value you want to change.
    :param lat: A 1D array of latitude value you want to change.
    :return:  The data after resolution processing has the same dimension as d.
    d包含的经纬度范围>=输入的lon,lat范围

        初始的经纬度与目标经纬度定义格式要相同，比如： 经度都是[-80,80]而不能一个是[-180，180],另一个是[0, 360]
    即：数组中各数值的位置可以不一一对应，但必须有一样的数值，比如西经40°，在两个经度数据中都用-40来表示

    '''
    # 先将原始数据的整体范围与目标范围缩小到一致
    max_lat_value = np.amax(lat)
    min_lat_value = np.amin(lat)
    max_lon_value = np.amax(lon)
    min_lon_value = np.amin(lon)
    max_lat0_index = list(lat0).index(max_lat_value)  # 获取原始经纬度与目标经纬度对应边界值的索引
    min_lat0_index = list(lat0).index(min_lat_value)
    max_lon0_index = list(lon0).index(max_lon_value)
    min_lon0_index = list(lon0).index(min_lon_value)
    a1 = min(max_lat0_index, min_lat0_index)  # 两个纬度可能排序不同，如[-80,80]或[80,-80]，选取较大索引
    a2 = max(max_lat0_index, min_lat0_index)
    b1 = min(max_lon0_index, min_lon0_index)
    b2 = max(max_lon0_index, min_lon0_index)
    data1 = d[:, a1:a2 + 1, b1:b2 + 1]  # 从d中截取与目标经纬度起始范围相同的数据集data1
    lat1 = lat0[a1:a2 + 1]
    lon1 = lon0[b1:b2 + 1]
    data = np.zeros((d.shape[0], len(lat), len(lon)))
    for i in range(len(lat)):
        y = lat[i]
        y1 = list(lat1).index(y)
        for j in range(len(lon)):
            x = lon[j]
            x1 = list(lon1).index(x)
            data[:, i, j] = data1[:, y1, x1]
    return data


u1_i = filter3D_lonlat(u_i, lon0, lat0, lon1, lat1)
u1_e = filter3D_lonlat(u_e, lon0, lat0, lon1, lat1)
u1_H = filter3D_lonlat(u_H, lon0, lat0, lon1, lat1)
v1_i = filter3D_lonlat(v_i, lon0, lat0, lon1, lat1)
v1_e = filter3D_lonlat(v_e, lon0, lat0, lon1, lat1)
v1_H = filter3D_lonlat(v_H, lon0, lat0, lon1, lat1)

lon1 = np.array(lon1)
lat1 = np.array(lat1)

levels1 = np.arange(-0.6, 0.66, 0.06)
levels2 = np.arange(-0.5, 0.55, 0.05)
levels3 = np.arange(-0.45, 0.5, 0.05)

fig1, ax1 = plt.subplots(3, 3, figsize=(11, 12))
contourf1 = ax1[0, 0].contourf(abs(lon), abs(lat), rho_i[0, :, :], levels1, cmap='coolwarm', extend='both')
ax1[0, 1].contourf(abs(lon), abs(lat), rho_e[0, :, :], levels1, cmap='coolwarm', extend='both')
ax1[0, 2].contourf(abs(lon), abs(lat), rho_H[0, :, :], levels1, cmap='coolwarm', extend='both')
contourf2 = ax1[1, 0].contourf(abs(lon), abs(lat), rho_i[l1, :, :], levels2, cmap='coolwarm', extend='both')
ax1[1, 1].contourf(abs(lon), abs(lat), rho_e[l1, :, :], levels2, cmap='coolwarm', extend='both')
ax1[1, 2].contourf(abs(lon), abs(lat), rho_H[l1, :, :], levels2, cmap='coolwarm', extend='both')
contourf3 = ax1[2, 0].contourf(abs(lon), abs(lat), rho_i[l2, :, :], levels2, cmap='coolwarm', extend='both')
ax1[2, 1].contourf(abs(lon), abs(lat), rho_e[l2, :, :], levels2, cmap='coolwarm', extend='both')
ax1[2, 2].contourf(abs(lon), abs(lat), rho_H[l2, :, :], levels2, cmap='coolwarm', extend='both')

ax1[0, 0].set_title('(a)    $ ρ^{isQG}_{z=0}$', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[0, 1].set_title('(b)    $ ρ^{eSQG}_{z=0}$', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[0, 2].set_title('(c)    $ ρ^{HYCOM}_{z=0}$', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[1, 0].set_title('(d)    $ ρ^{isQG}_{z=-300}$', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[1, 1].set_title('(e)    $ ρ^{eSQG}_{z=-300}$', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[1, 2].set_title('(f)    $ ρ^{HYCOM}_{z=-300}$', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[2, 0].set_title('(g)    $ ρ^{isQG}_{z=-700}$', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[2, 1].set_title('(h)    $ ρ^{eSQG}_{z=-700}$', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[2, 2].set_title('(i)    $ ρ^{HYCOM}_{z=-700}$', font={'family': 'Arial', 'size': 12})  # 添加标题

ax1[0, 0].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [0, 0]窗口
ax1[0, 0].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位
ax1[0, 1].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [0, 1]
ax1[0, 1].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位
ax1[1, 0].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [1, 0]
ax1[1, 0].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位
ax1[1, 1].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [1, 1]
ax1[1, 1].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位
ax1[0, 2].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [0, 2]
ax1[0, 2].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位
ax1[1, 2].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [1, 2]
ax1[1, 2].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位
ax1[2, 0].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [1, 1]
ax1[2, 0].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位
ax1[2, 1].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [0, 2]
ax1[2, 1].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位
ax1[2, 2].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位  [1, 2]
ax1[2, 2].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位

ax1[0, 0].set_yticks([31, 34, 37, 40])  # 设置y轴刻度标签
ax1[0, 1].set_yticks([31, 34, 37, 40])
ax1[0, 2].set_yticks([31, 34, 37, 40])
ax1[1, 0].set_yticks([31, 34, 37, 40])
ax1[1, 1].set_yticks([31, 34, 37, 40])
ax1[1, 2].set_yticks([31, 34, 37, 40])
ax1[2, 1].set_yticks([31, 34, 37, 40])
ax1[2, 0].set_yticks([31, 34, 37, 40])
ax1[2, 2].set_yticks([31, 34, 37, 40])

cax1 = fig1.add_axes([0.92, 0.67, 0.015, 0.2])  # 创建颜色条的位置和大小
bar1 = plt.colorbar(contourf1, cax=cax1, orientation='vertical')
bar1.set_ticks([-0.5, -0.3, -0.1, 0.1, 0.3, 0.5])  # 设置颜色条刻度

cax2 = fig1.add_axes([0.92, 0.39, 0.015, 0.2])
bar2 = plt.colorbar(contourf2, cax=cax2, orientation='vertical')
bar2.set_ticks([-0.4, -0.2, 0, 0.2, 0.4])

cax3 = fig1.add_axes([0.92, 0.11, 0.015, 0.2])
bar3 = plt.colorbar(contourf3, cax=cax3, orientation='vertical')
bar3.set_ticks([-0.4, -0.2, 0, 0.2, 0.4])


# 画速度场矢量图（带箭头）
style = dict(width=0.005, headwidth=2.2, headlength=1, headaxislength=2, color='dimgrey')
ax1[0, 0].quiver(abs(lon1), abs(lat1), u1_i[0, :, :], v1_i[0, :, :], **style)
ax1[0, 2].quiver(abs(lon1), abs(lat1), u1_H[0, :, :], v1_H[0, :, :], **style)
ax1[0, 1].quiver(abs(lon1), abs(lat1), u1_e[0, :, :], v1_e[0, :, :], **style)
ax1[1, 0].quiver(abs(lon1), abs(lat1), u1_i[l1, :, :], v1_i[l1, :, :], **style)
ax1[1, 1].quiver(abs(lon1), abs(lat1), u1_e[l1, :, :], v1_e[l1, :, :], **style)
ax1[1, 2].quiver(abs(lon1), abs(lat1), u1_H[l1, :, :], v1_H[l1, :, :], **style)
ax1[2, 0].quiver(abs(lon1), abs(lat1), u1_i[l2, :, :], v1_i[l2, :, :], **style)
ax1[2, 1].quiver(abs(lon1), abs(lat1), u1_e[l2, :, :], v1_e[l2, :, :], **style)
ax1[2, 2].quiver(abs(lon1), abs(lat1), u1_H[l2, :, :], v1_H[l2, :, :], **style)

ax1[0, 0].invert_xaxis()  # 反转坐标轴
ax1[0, 1].invert_xaxis()
ax1[0, 2].invert_xaxis()
ax1[1, 0].invert_xaxis()
ax1[1, 1].invert_xaxis()
ax1[1, 2].invert_xaxis()
ax1[2, 1].invert_xaxis()
ax1[2, 0].invert_xaxis()
ax1[2, 2].invert_xaxis()

# 添加整个窗口的 x 轴和 y 轴标签
fig1.text(0.5, 0.04, 'Longitude', ha='center', font={'family': 'Arial', 'size': 16})
fig1.text(0.04, 0.5, 'Latitude', va='center', rotation='vertical', font={'family': 'Arial', 'size': 16})

plt.subplots_adjust(wspace=0.3, hspace=0.35)

# ************************** profile
rho_ii = d1['rhoi'][:]
rho_is = d1['rhos'][:]

rms_ii = np.zeros(len(z))
rms_is, rms_i, rms_e, rms_H = np.zeros(len(z)), np.zeros(len(z)), np.zeros(len(z)), np.zeros(len(z))
index = 88
def rms(data):
    df = data.flatten()
    n = len(df)
    return np.sqrt(np.sum(df**2)/n)

for k in range(len(z)):
    rms_ii[k] = rms(rho_ii[k, :, :])
    rms_is[k] = rms(rho_is[k, :, :])
    rms_i[k] = rms(rho_i[k, :, :])
    rms_e[k] = rms(rho_e[k, :, :])
    rms_H[k] = rms(rho_H[k, :, :])

levels1 = np.arange(-0.4, 0.44, 0.04)
levels2 = np.arange(-0.2, 0.2, 0.01)
fig1, ax1 = plt.subplots(3, 2, figsize=(8, 9))
contourf11 = ax1[0, 0].contourf(abs(lon), -z, rho_is[:, index, :], levels1, cmap='coolwarm', extend='both')
contourf12 = ax1[0, 1].contourf(abs(lon), -z, rho_ii[:, index, :], levels1, cmap='coolwarm', extend='both')
contourf13 = ax1[1, 0].contourf(abs(lon), -z, rho_i[:, index, :], levels1, cmap='coolwarm', extend='both')
contourf14 = ax1[1, 1].contourf(abs(lon), -z, rho_e[:, index, :], levels1, cmap='coolwarm', extend='both')
contourf15 = ax1[2, 0].contourf(abs(lon), -z, rho_H[:, index, :], levels1, cmap='coolwarm', extend='both')

ax1[0, 0].set_title('(a)  $ \\rho^{s}$', font={'family': 'Arial', 'size': 12})  # 添加标题
ax1[0, 1].set_title('(b)  $ \\rho^{i}$', font={'family': 'Arial', 'size': 12})
ax1[1, 0].set_title('(c)  $ \\rho^{isQG}$', font={'family': 'Arial', 'size': 12})
ax1[1, 1].set_title('(d)  $ \\rho^{eSQG}$', font={'family': 'Arial', 'size': 12})
ax1[2, 0].set_title('(e)  $ \\rho^{HYCOM}$', font={'family': 'Arial', 'size': 12})

ax1[0, 0].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位
ax1[0, 1].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))
ax1[1, 0].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))
ax1[1, 1].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))
ax1[2, 0].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))

ax1[0, 0].invert_xaxis()  # 反转坐标轴
ax1[0, 1].invert_xaxis()
ax1[1, 0].invert_xaxis()
ax1[1, 1].invert_xaxis()
ax1[2, 0].invert_xaxis()

ax1[0, 0].set_yticks([0, -500, -1000, -1500, -2000, -2500])  # 设置y轴刻度标签
ax1[0, 1].set_yticks([0, -500, -1000, -1500, -2000, -2500])
ax1[1, 0].set_yticks([0, -500, -1000, -1500, -2000, -2500])
ax1[1, 1].set_yticks([0, -500, -1000, -1500, -2000, -2500])
ax1[2, 0].set_yticks([0, -500, -1000, -1500, -2000, -2500])
ax1[2, 1].set_yticks([0, -500, -1000, -1500, -2000, -2500])

ax1[0, 0].tick_params(axis='both', which='both', length=0, labelsize=9.8)  # 去掉轴刻度线
ax1[0, 1].tick_params(axis='both', which='both', length=0, labelsize=9.8)
ax1[1, 0].tick_params(axis='both', which='both', length=0, labelsize=9.8)
ax1[1, 1].tick_params(axis='both', which='both', length=0, labelsize=9.8)
ax1[2, 0].tick_params(axis='both', which='both', length=0, labelsize=9.8)
ax1[2, 1].tick_params(axis='both', which='both', length=0, labelsize=9.8)

ax1[2, 1].plot(rms_is, -z, color='slategrey', linestyle=':', label='s')  #
ax1[2, 1].plot(rms_ii, -z, color='teal', linestyle='-.', label='i')  # 画垂直分层
ax1[2, 1].plot(rms_i, -z, color='lightcoral', label='isQG')  # 画垂直分层
ax1[2, 1].plot(rms_e, -z, color='darkseagreen', label='eSQG')  #
ax1[2, 1].plot(rms_H, -z, color='royalblue', label='HYCOM')  # 画垂直分层

ax1[2, 1].legend(loc='best')  # 添加图例

ax1[2, 1].set_title("(f)   Root Mean Square", font={'family': 'Arial', 'size': 12})  # 添加轴标签


cax = fig1.add_axes([0.92, 0.3, 0.02, 0.5])  # 创建颜色条的位置和大小
bar = plt.colorbar(contourf15, cax=cax, orientation='vertical')  # vertical
bar.set_ticks([-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3])  # 设置颜色条刻度

fig1.text(0.5, 0.04, 'Longitude', ha='center', font={'family': 'Arial', 'size': 16})
fig1.text(0.01, 0.5, 'Depth(m)', va='center', rotation='vertical', font={'family': 'Arial', 'size': 16})

plt.subplots_adjust(wspace=0.3, hspace=0.3)

# ********** correlation
x = 98  # 52.16°W
y = 88  # 33.52°N

Aui = np.corrcoef(u_i[:, y, :], u_H[:, y, :])
Aue = np.corrcoef(u_e[:, y, :], u_H[:, y, :])

cor_ui = np.zeros(len(z))
cor_ue = np.zeros(len(z))
for k in range(len(z)):
    cor_ui[k] = Aui[k, 38+k]
    cor_ue[k] = Aue[k, 38+k]

Avi = np.corrcoef(v_i[:, y, :], v_H[:, y, :])
Ave = np.corrcoef(v_e[:, y, :], v_H[:, y, :])

cor_vi = np.zeros(len(z))
cor_ve = np.zeros(len(z))
for k in range(len(z)):
    cor_vi[k] = Avi[k, 38+k]
    cor_ve[k] = Ave[k, 38+k]

Ari = np.corrcoef(rho_i[:, y, :], rho_H[:, y, :])
Are = np.corrcoef(rho_e[:, y, :], rho_H[:, y, :])

cor_ri = np.zeros(len(z))
cor_re = np.zeros(len(z))
for k in range(len(z)):
    cor_ri[k] = Ari[k, 38+k]
    cor_re[k] = Are[k, 38+k]

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

RMSE_ri = np.zeros(len(z))
RMSE_re = np.zeros(len(z))
RMSE_ui = np.zeros(len(z))
RMSE_ue = np.zeros(len(z))
RMSE_vi = np.zeros(len(z))
RMSE_ve = np.zeros(len(z))

for k in range(len(z)):
    RMSE_ri[k] = get_rmse(rho_H[k, y, :], rho_i[k, y, :])
    RMSE_re[k] = get_rmse(rho_H[k, y, :], rho_e[k, y, :])
    RMSE_ui[k] = get_rmse(u_H[k, y, :], u_i[k, y, :])
    RMSE_ue[k] = get_rmse(u_H[k, y, :], u_e[k, y, :])
    RMSE_vi[k] = get_rmse(v_H[k, y, :], v_i[k, y, :])
    RMSE_ve[k] = get_rmse(v_H[k, y, :], v_e[k, y, :])

fig1, ax = plt.subplots(1, 3, figsize=(11,5))
ax[0].plot(cor_ri[0:33], -z[0:33], color='orchid', label='isQG')  # 画垂直分层, lightcoral
ax[0].plot(cor_re[0:33], -z[0:33], color='cornflowerblue', label='eSQG')  # royalblue
ax[0].legend(loc='best')  # 添加图例
ax[0].set_title(" (a)  Density", font={'family': 'Arial', 'size': 16})  # 添加轴标签


ax[1].plot(cor_ui[0:33], -z[0:33], color='lightsalmon', label='isQG')  # 画垂直分层, lightcoral
ax[1].plot(cor_ue[0:33], -z[0:33], color='slateblue', label='eSQG')  # rosybrown
ax[1].legend(loc='best')  # 添加图例
ax[1].set_title(" (b)  Zonal velocity ", font={'family': 'Arial', 'size': 16})  # 添加轴标签

ax[2].plot(cor_vi[0:33], -z[0:33], color='darkcyan', label='isQG')  # 画垂直分层, lightcoral, sandybrown
ax[2].plot(cor_ve[0:33], -z[0:33], color='darkblue', label='eSQG')  # slategray
ax[2].legend(loc='best')  # 添加图例
ax[2].set_title(" (c)  Meridional velocity ", font={'family': 'Arial', 'size': 16})  # 添加轴标签

fig1.text(0.5, 0.015, 'Correlation', ha='center', font={'family': 'Arial', 'size': 16})
fig1.text(0.03, 0.5, 'Depth(m)', va='center', rotation='vertical', font={'family': 'Arial', 'size': 16})
plt.subplots_adjust(wspace=0.3)
plt.show()
