import numpy as np
import sympy as sp


class LinearSearch:
    '''
    线性搜索方法的整合
    '''
    def gold(self,a,b,epsilon = 1e-3,f = lambda x:x):
        '''
        直接搜索---0.618法
        :param a: 区间左端点
        :param b: 区间右端点
        :param epsilon: 精度阈值
        :param f: 函数表达式,可以是lambda表达式
        :return: 收敛到精度的最优解
        '''
        t = 0.618
        al,ar = a + (1 - t) * (b - a), a + t * (b - a)
        fl,fr = f(al), f(ar)
        while abs(b - a) >= epsilon:
            if fl < fr:
                b = ar
                ar = al
                fr = fl
                al = (b - a) * (1 - t) + a
                fl = f(al)
            else:
                a = al
                al = ar
                fl = fr
                ar = (b - a) * t + a
                fr = f(ar)
        if fl < fr:
            return al
        else:
            return ar

    def stepinback(self,a0,y = 2,h = 0.1,f = lambda x: x,returnC = False):
        '''
        直接搜索---进退法
        :param a0: 给定初始点
        :param y: 步长前面乘的
        :param h: 步长
        :param f: 函数表达式
        :return: (returnC == True)a,c,b f(a)>f(c)<f(b)\
        (returnC == False) a,b
        '''
        a,b,c = None,None,None
        a1 = a0 + h
        if f(a0) > f(a1):
            a2 = a1 + y * h
            if f(a2) > f(a1):
                a,c,b = a0,a1,a2
            else:
                k = 2
                while f(a2) <= f(a1):
                    a0 = a1
                    a1 = a2
                    a2 = a2 + (y ** k) * h
                    k +=1
                a,c,b = a0,a1,a2
        else:
            t = a0
            a0 = a1
            a1 = t
            a2 = a1 - y * h
            if f(a2) > f(a1):
                a,c,b = a2,a1,a0
            else:
                k = 2
                while f(a2) <= f(a1):
                    a0 = a1
                    a1 = a2
                    a2 = a2 - (y ** k) * h
                    k += 1
                a,c,b = a2, a1, a0

        if returnC:
            return a,c,b
        else:
            return a,b

    def Goldstein(self,u,beta1,beta2,f = lambda x:x):
        '''
        Goldstein搜索法
        :param u: 试探点
        :param beta1: 精度要求1
        :param beta2: 精度要求2
        :param f: 函数表达式(lambda表达式)
        :return: 最优解
        '''
        max_iters = int('1' + '0'*6)
        assert beta1>0 and beta1<1 and beta2>0 and beta2<1 and beta1 < beta2,r"精度设置错误，必须在0到1之间，且beta1<beta2"
        umin,umax = 0, np.inf
        # 将Lambda表达式转换为符号表达式
        x = sp.Symbol('x')
        f = sp.sympify(f(x))
        f_diff = sp.diff(f,x)

        for i in range(max_iters):
            if f.subs(x,u) > f.subs(x,0) + beta1 * f_diff.subs(x,0) * u:
                umax = u
            elif f.subs(x,u) >= f.subs(x,0) + beta2 * f_diff.subs(x,0) * u:
                break
            else:
                umin = u
            if umax < np.inf:
                u = (umin + umax) / 2
            else:
                u = 2 * umin

        return u

    def preciseLinearSearh(self,f:lambda x:x):
        # 将Lambda表达式转换为符号表达式
        x = sp.Symbol('x')
        f = sp.sympify(f(x))
        f_diff = sp.diff(f, x)
        f_diff2 = sp.diff(f_diff,x)
        solutions = np.array(sp.solveset(f_diff,x))
        if solutions.shape[0] <= 0:
            raise "No Solution"
        elif solutions.shape[0] > 1:
            for x0 in solutions:
                if f_diff2.subs(x,x0) > 0:
                    return x0
        else:
            return solutions

class DescendNC:
    '''
    下降算法整合
    '''
    def __init__(self,kernal = "0.618" ):
        '''
        :param kernal:确定下降步长的算法,默认方法0.618法
        '''
        assert kernal in ('0.618','precise','goldstein'),r'目前计算包只支持0.618法(0.618),精确线性搜索(precise),goldstein法(goldstein)'
        self.kernal = kernal

    def speedDescend(self,f = lambda x,y: x+y,*symbols):

        grad = np.zeros(len(symbols))
        for symbol in symbols:
            pass












