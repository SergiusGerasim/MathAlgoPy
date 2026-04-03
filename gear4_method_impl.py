import numpy as np
import plotly.graph_objects as go
from tqdm import trange

def graph(data_accuracy, name = '', flag = True):
    if flag:
        layout_accuracy = go.Layout()
        fig1 = go.Figure(data=data_accuracy, layout=layout_accuracy)
        fig1.update_layout(yaxis_title="f(x)", xaxis_title="x", title=name)
        fig1.show()

def Gear(h, a, N, dydx, SOL,lamb, function_count):
    Y = np.zeros((function_count, N))
    x = np.zeros(N)
    
    x[0] = a
    Y[:,0] = SOL


    for i in range(1, 4):

        k = np.array([np.zeros(4) for _ in range(function_count)])

        k[:,0] = dydx(x[i-1], *[Y[n][i-1] for n in range(function_count)])

        for j in range(1, 4):
            for n in range(function_count):
                if j == 3:
                    k[n][j] = dydx(x[i-1] + h, *[Y[m][i-1] + h*k[m][j-1] for m in range(function_count)])[n]
                else:
                    k[n][j] = dydx(x[i-1] + 0.5*h, *[Y[m][i-1] + 0.5*h*k[m][j-1] for m in range(function_count)])[n]
        
        x[i] = x[i-1] + h
        
        Y[:,i] = Y[:,i-1] + (h / 6) * (k[:,0] + 2*k[:,1] + 2*k[:,2] + k[:,3])
        


    for i in trange(4, N):
        x[i] = x[i-1] + h

        Y[:,i] = ( 48*Y[:,i-1] - 36*Y[:,i-2] + 16*Y[:,i-3] - 3*Y[:,i-4] ) / (
            25 - 12 * lamb * h
        )

    return x, Y

def Gear_iter(h, a, N, dydx, SOL, function_count):
    Y = np.zeros((function_count, N))
    x = np.zeros(N)
    
    x[0] = a
    Y[:,0] = SOL


    for i in range(1, 4):

        k = np.array([np.zeros(4) for _ in range(function_count)])

        k[:,0] = dydx(x[i-1], *[Y[n][i-1] for n in range(function_count)])

        for j in range(1, 4):
            for n in range(function_count):
                if j == 3:
                    k[n][j] = dydx(x[i-1] + h, *[Y[m][i-1] + h*k[m][j-1] for m in range(function_count)])[n]
                else:
                    k[n][j] = dydx(x[i-1] + 0.5*h, *[Y[m][i-1] + 0.5*h*k[m][j-1] for m in range(function_count)])[n]
        
        x[i] = x[i-1] + h
        
        Y[:,i] = Y[:,i-1] + (h / 6) * (k[:,0] + 2*k[:,1] + 2*k[:,2] + k[:,3])
        


    for i in trange(4, N):
        x[i] = x[i-1] + h
        
        y0 = (1/25) * ( 12*h*dydx(x[i-1], *[Y[n][i-1] for n in range(function_count)]) 
                       + 48*Y[:,i-1] - 36*Y[:,i-2] + 16*Y[:,i-3] - 3*Y[:,i-4] ) 
        y1 = (1/25) * ( 12*h*dydx(x[i], *y0) 
                       + 48*Y[:,i-1] - 36*Y[:,i-2] + 16*Y[:,i-3] - 3*Y[:,i-4] ) 
 
        yn = np.array([y0, y1])
        it = 0
        while(np.linalg.norm(y1 - y0, ord=np.inf)/np.linalg.norm(y0,ord=np.inf) > 10**(-5)):
            it+=1
            y0 = y1
            y1 = (1/25) * ( 12*h*dydx(x[i], *y0 ) 
                       + 48*Y[:,i-1] - 36*Y[:,i-2] + 16*Y[:,i-3] - 3*Y[:,i-4] )
        print(f"шаг {i}, кол-во итераций {it}")
        Y[:,i] = y1

    return x, Y