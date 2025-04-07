# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 15:28:20 2013

@author: akels
"""
from __future__ import division, print_function
from os import sys
sys.path.append('cpresources')
from pylab import *


h = 2e-18*10
hbar = 1.0546e-36
L = 1e-8
M = 9.109e-31
N = 1000 # Grid slices

a = L/N

def complex_arg(trans):
    """装饰器函数，用于处理复数参数的变换"""
    def f(y):
        # 分别对实部和虚部进行变换
        return trans(real(y)) + 1j*trans(imag(y))
    return f

@complex_arg
def dst(y):
    """离散正弦变换(DST)
    
    参数:
        y (array): 输入数组，可以是实数或复数
        
    返回:
        array: 变换后的系数数组
    """
    N = len(y)
    # 准备扩展数组用于FFT计算
    y2 = empty(2*N,float)
    y2[0] = y2[N] = 0.0  # 设置边界条件
    y2[1:N] = y[1:]      # 填充前半部分
    y2[:N:-1] = -y[1:]   # 填充后半部分(奇对称)
    
    # 计算FFT并取虚部作为DST结果
    a = -imag(rfft(y2))[:N]
    a[0] = 0.0  # 确保DC分量为0
    
    return a


######################################################################
# 1D inverse DST Type-I

@complex_arg
def idst(a):
    """一维离散正弦逆变换(DST-I)
    
    参数:
        a (array): 输入数组，包含DST系数
        
    返回:
        array: 逆变换后的实数/复数数组
    """
    N = len(a)
    # 准备复数数组用于逆FFT
    c = empty(N+1,complex)
    c[0] = c[N] = 0.0  # 边界条件设置为0
    c[1:N] = -1j*a[1:]  # 填充中间系数
    
    # 执行逆实数FFT并截取前N个点
    y = irfft(c)[:N]
    y[0] = 0.0  # 确保边界条件

    return y



ksi = zeros(N+1,complex)  # 初始化波函数数组

def ksi0(x):
    """初始波函数"""
    x0 = L/2       # 波包中心位置
    sigma = 1e-10  # 波包宽度
    k = 5e10       # 波数
    # 高斯波包乘以平面波
    return exp(-(x-x0)**2/2/sigma**2)*exp(1j*k*x)

# 初始化空间网格和波函数
x = linspace(0,L,N+1)  # 空间坐标
ksi[:] = ksi0(x)       # 设置初始波函数
ksi[[0,N]]=0           # 边界条件设为0

# 计算初始时刻的DST系数
b0 = dst(ksi)

t = 1e-18
b_ = b0*exp(1j*pi**2*hbar*arange(1,N+2)**2/2/M/L**2*t)

ksi_ = idst(b_)
plot(ksi_)
show()


# 可视化部分
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# 创建3D图形
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 初始化曲线数据
line, = ax.plot(x-L/2, real(ksi)*1e-9, imag(ksi)*1e-9)

# 设置坐标轴标签
ax.set_xlabel('Position')
ax.set_ylabel('Real part')
ax.set_zlabel('Imaginary part')

# 主循环 - 时间演化
t = 0
for _ in range(1000):  # 限制循环次数
    plt.pause(0.03)  # 控制动画帧率
    
    # 计算时间演化后的系数
    b_ = b0*exp(1j*pi**2*hbar*arange(1,N+2)**2/2/M/L**2*t)
    
    # 逆变换得到实时波函数
    ksi_ = idst(b_)
    
    # 更新曲线数据
    line.set_data(x-L/2, real(ksi_)*1e-9)
    line.set_3d_properties(imag(ksi_)*1e-9)
    
    t += h*5  # 时间步进

plt.show()