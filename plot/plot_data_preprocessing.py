import numpy as np
import netCDF4 as nc
import matplotlib.ticker as mticker
from numpy import pi, sin, zeros, arange
from pylab import meshgrid

import matplotlib.pyplot as plt


def anomaly(lon, lat, var):
    """
    Compute the anomaly of a variable with respect to its spatial mean.

    Parameters
    ----------
    lon : Array of longitudes.
    lat : Array of latitudes.
    var : Array of variable values.

    Returns
    -------
    anomaly : array-like
        Anomaly of the variable with respect to its spatial mean.
    """

    if lon.ndim == 1:
        y1, x1 = meshgrid(lon, lat)
    else:
        y1, x1 = lon, lat

    if var.ndim == 2:
        coef = fit2Dsurf(x1, y1, var)[0]
        varm = fit2poly(coef, x1, y1)
        return var - varm

    elif var.ndim == 3:
        tmp = zeros(var.shape)
        for i in arange(var.shape[0]):
            coef = fit2Dsurf(x1, y1, var[i, :, :])[0]  # 拟合系数
            varm = fit2poly(coef, x1, y1)  # 拟合函数
            tmp[i, :, :] = var[i, :, :] - varm
        return tmp


def fit2Dsurf(x, y, p):
    """
    Given y0=f(t0), find the best fit
    p = a + bx + cy + dx**2 + ey**2 + fxy
    and return a, b, c, d, e, f

    :param x: 1D array of x values
    :param y: 1D array of y values
    :param p: 1D array of z values
    :return: a tuple of coefficients (a, b, c, d, e, f)
    """
    from scipy.optimize import leastsq

    def fit2poly(params, x, y, p):
        a, b, c, d, e, f = params
        err = p - (a + b * x + c * y + d * x**2 + e * y**2 + f * x * y)
        return err

    x = x.flatten()
    y = y.flatten()
    p = p.flatten()
    params = np.array([p.mean(), 1e-3, 1e-3, 1e-6, 1e-6, 1e-6])
    coefs = leastsq(fit2poly, params, args=(x, y, p))
    return coefs


def fit2poly(params, x, y):
    """
    Returns the polynomial fit p = a + bx + cy + dx**2 + ey**2 + fxy
    given the coefficients c and input data x and y.

    Parameters:
        params (list): List of polynomial coefficients
        x (ndarray): Array of x-values
        y (ndarray): Array of y-values

    Returns:
        fit (ndarray): Array of polynomial fit values
    """
    a, b, c, d, e, f = params
    fit = (a + b*x + c*y + d*x**2 + e*y**2 + f*x*y)
    return fit


d = nc.Dataset('path_GS.nc')
lat = d['lat'][:]  # 26N-40N, 351
lon = d['lon'][:]
z = d['depth'][:]  # 0-1000m, 33
N2 = d['N2'][:]
ssd = d['ssd'][:]  # 33 * 176 * 351
ssdl = d['ssdl'][:]
ssda = d['ssda'][:]
M = round(np.mean(ssd), 1)
print(M)
ssh = d['ssh'][:]  # 176 * 351
sshl = d['sshl'][:]
ssha = d['ssha'][:]

levels1= np.arange(-0.8, 0.9, 0.1)
levels2 = np.arange(-0.5, 0.6, 0.1)
levels3= np.arange(-0.6, 0.7, 0.1)
"""  画图———密度异常与高度异常（根据实验前的数据画图）  """

fig1, ax = plt.subplots(1, 3, figsize=(12, 4))
contourf1 = ax[0].contourf(abs(lon), abs(lat), ssd-M, levels1, cmap='coolwarm', extend='both')  # ssd anomaly
contour1 = ax[0].contour(abs(lon), abs(lat), ssh, inline=True, fontsize=6, colors='w', linewidths=0.6)  # ssh anomaly
ax[0].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位
ax[0].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位
ax[0].set_title('(a)', font={'family': 'Arial', 'size': 16})  # 添加标题
ax[0].set_xlabel('+1026.2 kg/$m^3$', labelpad=-2, x=0.9, font={'family': 'Arial', 'size': 9})  # 添加标题
ax[0].invert_xaxis()  # 反转坐标轴
ax[0].set_yticks([31, 34, 37, 40])

contourf2 = ax[1].contourf(abs(lon), abs(lat), ssdl-M, levels2, cmap='coolwarm', extend='both')  # ssd anomaly
contour2 = ax[1].contour(abs(lon), abs(lat), sshl, inline=True, fontsize=6, colors='w', linewidths=0.6)  # ssh anomaly
ax[1].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位
ax[1].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位
ax[1].set_title('(b)', font={'family': 'Arial', 'size': 16})  # 添加标题
ax[1].set_xlabel('+1026.2 kg/$m^3$', labelpad=-2, x=0.9, font={'family': 'Arial', 'size': 9})  # 添加标题
ax[1].invert_xaxis()  # 反转坐标轴
ax[1].set_yticks([31, 34, 37, 40])

contourf3 = ax[2].contourf(abs(lon), abs(lat), ssda, levels3, cmap='coolwarm', extend='both')  # ssd anomaly
contour3 = ax[2].contour(abs(lon), abs(lat), ssha, inline=True, fontsize=6, colors='w', linewidths=0.6)  # ssh anomaly
ax[2].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°W'))  # 给经度加单位
ax[2].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.0f°N'))  # 给纬度加单位
ax[2].set_title('(c)', font={'family': 'Arial', 'size': 16})  # 添加标题
ax[2].invert_xaxis()  # 反转坐标轴
ax[2].set_yticks([31, 34, 37, 40])

bar1 = plt.colorbar(contourf1, orientation='horizontal', pad=0.13, shrink=0.9)
bar1.set_ticks([-0.7, -0.3, 0, 0.4])  # 给colorbar添加刻度标签
bar2 = plt.colorbar(contourf2, orientation='horizontal', pad=0.13, shrink=0.9)
bar2.set_ticks([-0.3, 0, 0.4])
bar3 = plt.colorbar(contourf3, orientation='horizontal', pad=0.13, shrink=0.9)
bar3.set_ticks([-0.3, 0, 0.4])

fig1.text(0.5, 0.05, 'Longitude', ha='center', font={'family': 'Arial', 'size': 16})
fig1.text(0.05, 0.5, 'Latitude', va='center', rotation='vertical', font={'family': 'Arial', 'size': 16})
plt.subplots_adjust(wspace=0.3)  # hspace=0.5 上下间距

"""  画图———浮力频率分层剖面  """
d1 = nc.Dataset('path.nc')
n = d1['N2'][:]  # 表层=1-19的均值
N2 = d1['N20'][:]
lat = d1['lat'][:]  # 26N-40N, 351
f0 = 4.0*pi/(24*3600)*sin(lat.mean()*pi/180.0)
N = np.sqrt(N2)
nb = N/abs(f0)  # prandtl ratio
N0 = np.mean(N[:])  # 常数N
n0 = N0/abs(f0)
A = np.zeros(len(z))
B = n0+A
fig1, ax = plt.subplots(1, 2, figsize=(10, 6))
ax[0].plot(N2*10**5, -z)  # 画垂直分层
ax[0].plot(n[0:13] * 10 ** 5, -z[0:13], color='slategrey', linestyle='--')
# ax1.set_xlim(3, 6.5)  # Y轴范围2-7
ax[0].set_title("(a) $N^{2}(z)$", font={'family': 'Arial', 'size': 16})  # 添加轴标签
ax[0].set_xlabel('($\\times  10^{-5}$)', labelpad=-15, x=1.04, font={'family': 'Arial', 'size': 9})  # 添加标题
ax[1].plot(nb, -z, color='dimgrey', linestyle='-.')  # 画垂直分层
ax[1].plot(B, -z)
ax[1].set_title("(b) $N_0/|f_0|$", font={'family': 'Arial', 'size': 16})  # 添加轴标签
plt.subplots_adjust(wspace=0.3)
fig1.text(0.02, 0.5, 'Depth(m)', va='center', rotation='vertical', font={'family': 'Arial', 'size': 16})

plt.show()
