# 【例題2.8】

## (6) 洪水追跡計算
### ①Kinematic wave法による洪水追跡計算プログラム

- kinematic wave 法による洪水追跡計算を行います．

- 計算条件は，エクセルファイル”KWInput.xlsx”をご確認ください．
    - 計算条件の設定は，"計算条件"のシートで行います．
    - 上流端流量の設定は，"上流端流量の経時変化"のシートで行います．
    - 上流端流量ハイドログラフは，プログラム"Discharge/discharge.py"で作成します．


- 各時刻の水深，流積，流速，流量の空間分布が，"KW/ResKWave.xlsx"に出力されます．

- "KW/28_8_6graph.py"で，計算結果をグラフにできます．


```python
'''kinematic wave 法による洪水追跡計算を行うプログラム'''
# !/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import scipy.optimize
import decimal
# 定数の設定
GRAVITY_ACCELERATION: float = 9.81  # 重力加速度
INPUT_FILE_NAME = "KWInput.xlsx"
SHEET1_NAME = "計算条件"
SHEET2_NAME = "上流端流量の経時変化"
OUTPUT_FILE_NAME = "KW/ResKWave.xlsx"
EPS_DEPTH = 1e-15
decimal.getcontext().prec = 10
```

#### 入力データを読み込むクラス


```python
class SetData:
    """入力データを読み込むクラス"""
    def __init__(self):
        """データの読み込み"""
        self.data = pd.read_excel(INPUT_FILE_NAME, sheet_name=None)
        # 計算条件の読み込み
        print(self.data[SHEET1_NAME])
        __colName = self.data[SHEET1_NAME].columns
        # 河床勾配Ib
        self.Ib = self.data[SHEET1_NAME][__colName[1]][0]
        # 川幅B(m)
        self.B = self.data[SHEET1_NAME][__colName[1]][1]
        # 粗度係数n
        self.n = self.data[SHEET1_NAME][__colName[1]][2]
        # 計算距離
        self.xEnd = self.data[SHEET1_NAME][__colName[1]][3]*1000.
        # 分割数
        self.totalDivNo = int(self.data[SHEET1_NAME][__colName[1]][4])
        # 全格子数
        self.totalGridNo = self.totalDivNo+1
        # 格子幅dx
        self.dx = self.xEnd/self.totalDivNo
        # 計算時間 hr ->sec
        self.timeEnd = self.data[SHEET1_NAME][__colName[1]][5]*3600.
        # 出力時間 min -> sec
        self.outputTime = self.data[SHEET1_NAME][__colName[1]][6]*60.0
        # クーラン数
        self.Cr = self.data[SHEET1_NAME][__colName[1]][7]
        # 上流端流量の読み込み
        print(self.data[SHEET2_NAME])
        __colName = self.data[SHEET2_NAME].columns
        self.tb = np.array(self.data[SHEET2_NAME][__colName[0]])*3600.  # 時刻
        self.Qb = np.array(self.data[SHEET2_NAME][__colName[1]])  # 流量
        # 水深h,流量Qの配列を作成
        self.h = np.zeros(self.totalDivNo+1)
        self.Q = np.zeros(self.totalDivNo+1)

    def setInitCondition(self):
        """初期条件の設定"""
        self.x = np.array([self.dx*i for i in range(self.totalGridNo)])
        self.zb = self.Ib*(self.xEnd-self.x)
        self.Q[:] = self.Qb[0]
        self.h[:] = FlowParam.calch0(self.n, self.B, self.Ib, self.Qb[0])
        self.A = self.B*self.h
        self.R = FlowParam.calcR(self.A, FlowParam.calc_s(self.h, self.B))
        self.U = self.Q/self.A
        self.H = self.h+self.zb
        # データフレームの作成
        self.df = pd.DataFrame(self.x/1000.0, columns=['x(km)'])
        self.df['zb(m)'] = self.zb
```

#### 流れの諸量を計算するクラス


```python
class FlowParam:
    """流れの諸量を計算するクラス"""
    # 流積
    @staticmethod
    def calcA(_h, _B):
        return(_h*_B)

    # 潤辺
    @staticmethod
    def calc_s(_h, _B):
        return(_B+2.*_h)

    # 径深
    @staticmethod
    def calcR(_A, _s):
        return(_A/_s)

    # 等流時の流量を求める関数
    @staticmethod
    def calcQ0(_n, _A, _R, _Ib):
        return(1./_n*_A*_R**(2./3.)*_Ib**.5)

    # 等流の関係
    @staticmethod
    def __eqUniform(_h, _n, _B, _Ib, _Q):
        __A = FlowParam.calcA(_h, _B)
        __s = FlowParam.calc_s(_h, _B)
        __R = FlowParam.calcR(__A, __s)
        return(np.abs(-_Ib+_n**2.*(_Q/__A)**2./__R**(4./3.)))

    # 等流水深を求める関数
    @staticmethod
    def calch0(_n, _B, _Ib, _Q):
        __h0 = scipy.optimize.fmin(FlowParam.__eqUniform,
                                   x0=[(_n*_n*(_Q/_B)**2./_Ib)**(3./10.)],
                                   xtol=EPS_DEPTH,  disp=False,
                                   args=(_n, _B, _Ib, _Q, ))[0]
        return(__h0)
```

#### Kinematic wave 法のクラス

- 基礎方程式

$\frac{\partial A}{\partial t}+\frac{\partial Q}{\partial t} = 0$

$Q = \frac{1}{n} A R^{2/3}I_b^{1/2}$

- 離散化式

$A_i^{t+\Delta t} = A_i^{t}-\frac{\Delta t}{\Delta x}\left(Q_i^{t}-Q_{i-1}^t\right) $

- 境界条件式

$I_b-\frac{n^2}{{R_0^{t+\Delta t}}^{4/3}}\left(\frac{Q_0^{t+\Delta t}}{A_0^{t+\Delta t}}\right)^2 =0 $

上式より， $A_0^{t+\Delta t}$ を算出．


```python
class KinematicWave(SetData):
    """kinematic wave 法のクラス"""
    def __init__(self):
        # 変数等の設定
        super().__init__()
        super().setInitCondition()

    # ファイルへの書き出し
    def __writeFile(self, _fileNo, _time):
        # 計算結果の整理
        self.df['H(m)'] = self.h+self.zb
        self.df['h(m)'] = self.h
        self.df['A(m2)'] = self.A
        self.df['Q(m3/s)'] = self.Q
        self.df['U(m/s)'] = self.Q/self.A
        # sheet名
        __time = decimal.Decimal(_time)/decimal.Decimal(3600.0)
        __name = str(__time)+" hr"
        # ファイルへの出力
        if(_fileNo == 0):
            with pd.ExcelWriter(OUTPUT_FILE_NAME, mode='w') as writer:
                self.df.to_excel(writer, sheet_name=__name, index=False)
        else:
            with pd.ExcelWriter(OUTPUT_FILE_NAME, mode='a') as writer:
                self.df.to_excel(writer, sheet_name=__name, index=False)

    # 時間の刻み幅dtの計算
    def __calc_dT(self, _outputTime, _lambda):
        __dT = np.min([self.Cr*self.dx/_lambda, self.outputTime-_outputTime])
        return(__dT)

    # lambdaの計算
    def __calcLambda(self, _A, _Q, _R_B):
        return(np.max((5./3.-4./3.*_R_B)*(_Q/_A)))

    # Aの更新
    def __calcNewA(self, _At, _Qt, _dT):
        return(_At[1:]-_dT/self.dx*(_Qt[1:]-_Qt[:-1]))

    # 上流端Qの更新
    def __calcNewBC(self, _time):
        __i = np.where(self.tb > _time)[0][0]
        __dQ_dT = (self.Qb[__i]-self.Qb[__i-1])/(self.tb[__i]-self.tb[__i-1])
        __Qb = self.Qb[__i-1]+__dQ_dT*(_time-self.tb[__i-1])
        __hb = FlowParam.calch0(self.n, self.B, self.Ib, __Qb)
        __Ab = FlowParam.calcA(__hb, self.B)
        return(__Qb, __hb, __Ab)

    # 水深，潤辺，径深，流量の更新
    def __updateValue(self):
        self.h = self.A/self.B
        self.s = FlowParam.calc_s(self.h, self.B)
        self.R = FlowParam.calcR(self.A, self.s)
        self.Q[1:] = FlowParam.calcQ0(self.n, self.A[1:], self.R[1:], self.Ib)

    def procKW(self):
        """計算手順"""
        __time = 0.0
        __outputTime = 0.0
        __fileNo = 0
        # 初期条件のファイルへの出力
        self.__writeFile(__fileNo, __time)
        # 計算開始
        while(1):
            # 時刻tのAとQ
            __At = np.copy(self.A)
            __Qt = np.copy(self.Q)
            __lambda = self.__calcLambda(__At, __Qt, self.R/self.B)
            # dTの計算
            __dT = self.__calc_dT(__outputTime, __lambda)
            # 時刻と出力時間の更新
            __time += __dT
            __outputTime += __dT
            # Aの更新
            self.A[1:] = self.__calcNewA(__At, __Qt, __dT)
            # 上流端の流量と流積の更新
            (self.Q[0], self.h[0], self.A[0]) = self.__calcNewBC(__time)
            # 水深，潤辺，径深，流量の更新
            self.__updateValue()
            # ファイルへの出力
            if (__outputTime >= self.outputTime):
                print(__time)
                __fileNo += 1
                self.__writeFile(__fileNo, __time)
                __outputTime = 0.0
            # 計算時間を超えたら終了
            if (__time >= self.timeEnd):
                break
```

#### main関数


```python
if __name__ == "__main__":
    """main関数"""
    kw = KinematicWave()  # データの読み込み
    kw.procKW()  # 計算の実行
    print("Simulation is Over !")
```

                  変数名        数値
    0          河床勾配Ib    0.0005
    1          川幅B(m)  200.0000
    2  粗度係数n(m-1/3・s)    0.0300
    3        計算距離(km)  100.0000
    4        計算区間の分割数  250.0000
    5        計算時間（hr）   48.0000
    6       出力時間(min)   20.0000
    7         クーラン数Cr    0.8500
             t(hr)  Q(m3/s)
    0     0.000000    200.0
    1     0.166667    200.0
    2     0.333333    200.0
    3     0.500000    200.0
    4     0.666667    200.0
    ..         ...      ...
    296  49.333333    200.0
    297  49.500000    200.0
    298  49.666667    200.0
    299  49.833333    200.0
    300  50.000000    200.0
    
