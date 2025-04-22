#IMPORTANTO BIBLIOTECAS
from scipy.integrate import solve_ivp      #resolucao de E.D.
import numpy as np                         #computação científica e análise de dados (operacoes matematicas)
import matplotlib.pyplot as plt            #plotar graficos
import math
from flask import Flask, render_template, request
import io
import base64

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('simulador.html')

@app.route('/simulador')
def simulador():
    return render_template('simulador.html')

@app.route('/teoria')
def teoria():
    return render_template('teoria.html')

@app.route('/ajuda')
def ajuda():
    return render_template('ajuda.html')

@app.route('/autores')
def autores():
    return render_template('autores.html')

@app.route('/simular', methods=['GET','POST'])
def simular():
    if request.method == 'POST':
      # Pegando dados do formulário
      tempo = float(request.form['tempo'])
      step = float(request.form['step'])
      diametro = float(request.form['diametro']) * 1000  # Converter para mm
      curso = float(request.form['altura'])
      cilindrada_total = float(request.form['volume'])
      num_cilindros = int(request.form['num_cilindros'])
      taxa_comp = float(request.form['taxa_comp'])
      Pressao_entrada = float(request.form['pressao_entrada'])
      Temperatura_entrada = float(request.form['temperatura_entrada'])
      rotacao_inicial = float(request.form['rotacao'])

      #TEMPO DE SIMULACAO

      metodo = 'BDF'  # RK45 RK23 DOP853 BDF

      #DEFINICAO DOS PARAMETROS DO MODELO

      Ueo=213040                   # energia interna especifica inicial J/kg
      Po=101325                      # Pressao inicial Pa
      So=6863.05                     # entropia inicial J/kg.K
      Cv=722.85                      # calor especifico a volume constante J/kg.K
      Cp=1006                        # calor especifico a pressao constante J/kg
      k = Cp/Cv
      Q=5500                         # taxa de calor W

      #PARAMETROS AUXILIARES
      R2 = 0.0821 #constante universal dos gases em atm.L/mol.K
      Pressao_inicial = 1 #pressao do ar de entrada em atm
      R = 8.314 # Constante universal dos gases Pa.m^3/(mol.K)
      Temperatura_inicial = To = 298 #em K
      Ueo=213040                   # energia interna especifica inicial J/kg
      So=6863.05                     # entropia inicial J/kg.K
      Cpp = 29.1  #calor especifico a pressao constante
      M = 29 #massa molar do ar

      P2= 1823850 #em Pa

      #CALCULO DOS PARAMETROS DO MODELO

      Vdu =  cilindrada_total/num_cilindros
      V2  =  Vdu/(taxa_comp-1)   #volume morto em cm^3
      V1 = Vdu + V2
      v = 2*curso*rotacao_inicial        # velocidadde em m/s
      mo = (Pressao_inicial*Vdu*M)/(R2*Temperatura_inicial) #massa inicial em kg
      me = ((Pressao_entrada*M)/(R2*Temperatura_entrada))*((v*3.14*diametro)/4) #vazao massica de entrada em kg/s
      Uef=(Ueo - R*(Temperatura_entrada - Temperatura_inicial) + R*Temperatura_entrada*math.log(Pressao_entrada/Pressao_inicial))/M
      Veo= (R2*Temperatura_inicial)/(Pressao_inicial*M) #volume específico inicial em m^3/kg
      Vef= (R2*Temperatura_entrada)/(Pressao_entrada*M) #volume especifico final em m^3/kg
      Se= So + ((Cpp*math.log(Temperatura_entrada/Temperatura_inicial)-R*math.log(Pressao_entrada/Po)))/M #entropia de entrada
      W_rev=(P2*V2-Pressao_entrada*V1)/((k-1)*mo*M*tempo)
      W_res=0.312*W_rev

      # MODELAGEM DO 1 TEMPO - ADMISSAO

      def ode_func(t, y):   # t: variavel tempo em segundos e y: vetor de estados [massa, temperatura, entropia, exergia]

        m1 = y[0] # MASSA
        T1 = y[1] # TEMPERATURA
        S1 = y[2] # ENTROPIA
        X1 = y[3] # EXERGIA

        dm_dt = me
        dT_dt = (me/m1)*(k*Temperatura_entrada-T1)-(Q/(m1*Cv))
        dS_dt= (me/m1)*(Se-S1)
        dX_dt = (me/m1)*((v*v)/2)-X1

        return [dm_dt, dT_dt, dS_dt, dX_dt]

      #Definindo as condicoes iniciais

      m10=0.0002240  
      T10=298
      S10=6863.05
      X10=0

      CI=[m10,T10,S10,X10]

      #Tempo de simulacao
      t_span = (0, tempo)

      #comando para chamar o metodo de integracao
      sol = solve_ivp(ode_func, t_span, CI, method= metodo, max_step = step)    # RK45 RK23 DOP853 BDF

      # EXTRAINDO OS RESULTADOS
      t1 = sol.t
      m1 = sol.y[0]
      T1 = sol.y[1]
      S1 = sol.y[2]
      X1 = sol.y[3]

      # MODELAGEM DO 2 TEMPO: COMPRESSAO

      def ode_func2(t, y):   # t: variavel tempo em segundos e y: vetor de estados [massa, temperatura, entropia, exergia]

        m2 = y[0] # MASSA
        T2 = y[1] # TEMPERATURA
        S2 = y[2] # ENTROPIA
        X2 = y[3] # EXERGIA

        dm_dt = 0
        dT_dt = (W_rev - W_res)/(m2*Cv + (2*m2*R)/M)
        dS_dt= 0
        dX_dt = W_res/m2

        return [dm_dt, dT_dt, dS_dt, dX_dt]

      # Condições iniciais da compressão (usando o final da admissão)
      m20 = m1[-1]
      T20 = T1[-1]
      S20 = S1[-1]
      X20 = X1[-1]

      CI2 = [m20, T20, S20, X20]

      # Tempo de simulação para compressão
      t_span2 = (tempo, 2 * tempo)

      # Resolver a compressão
      sol2 = solve_ivp(ode_func2, t_span2, CI2, method=metodo, max_step=step)

      # Extraindo resultados
      t2 = sol2.t
      m2 = sol2.y[0]
      T2 = sol2.y[1]
      S2 = sol2.y[2]
      X2 = sol2.y[3]

      # MODELAGEM DO 3 TEMPO: COMBUSTAO

      def ode_func3(t, y):   # t: variavel tempo em segundos e y: vetor de estados [massa, temperatura, entropia, exergia]

        m3 = y[0] # MASSA
        T3 = y[1] # TEMPERATURA
        S3 = y[2] # ENTROPIA
        X3 = y[3] # EXERGIA

        MF = m3 - 1.8048e-3
        CT = 2.137e12*math.exp(-15098/T3)

        dm_dt = 0.4118 - (2.829e-5)*CT*((34.52*MF)**0.25)*(MF/0.00187)
        dT_dt = ((1430.8 - 1.7*T3 - (261.6*(T3-473)+377.7*(T3-1131)+1611.17*(T3-1131)) + (44.24e6)*(2.829e-5)*CT*((34.523*MF)**0.25)*((535.1*MF)**1.5)*0.0000584))/1.8852
        dS_dt= -(S3/MF)*(0.4118 - (2.8929e-5)*CT*((34.52*MF)**0.25)*((MF/1.8688e-3)**1.5)) + (56717.5/MF)
        dX_dt = -(X3/MF)*(0.4118 - 2.8929e-5)*CT*((34.52*MF)**0.25)*((MF/1.8688e-3)**1.5) - ((1.7*T3 - 14308)/MF) - (1679610/MF) + (321987.03e3/MF)

        return [dm_dt, dT_dt, dS_dt, dX_dt]

      #Definindo as condicoes iniciais

      m30=m2[-1]; 
      T30=T2[-1]; 
      S30=S2[-1]; 
      X30=X2[-1]

      CI3=[m30,T30,S30,X30]

      #Tempo de simulacao
      t_span3 = (2*tempo, 3*tempo)

      #comando para chamar o metodo de integracao
      sol3 = solve_ivp(ode_func3, t_span3, CI3, method=metodo, max_step = step)    # RK45 RK23 DOP853 BDF

      # EXTRAINDO OS RESULTADOS
      t3 = sol3.t
      m3 = sol3.y[0]
      T3 = sol3.y[1]
      S3 = sol3.y[2]
      X3 = sol3.y[3]

      # MODELAGEM DO 4 TEMPO: ESCAPE

      def ode_func4(t, y):   # t: variavel tempo em segundos e y: vetor de estados [massa, temperatura, entropia, exergia]

        m4 = y[0] # MASSA
        T4 = y[1] # TEMPERATURA
        S4 = y[2] # ENTROPIA
        X4 = y[3] # EXERGIA

        dm_dt = 0
        dT_dt = -10.4*T4 - 3.964
        dS_dt= 0
        dX_dt = 14226.64/m4

        return [dm_dt, dT_dt, dS_dt, dX_dt]

      #Definindo as condicoes iniciais

      m40=m3[-1]
      T40=T3[-1]
      S40=S3[-1]
      X40=X3[-1]

      CI4=[m40,T40,S40,X40]

      #Tempo de simulacao
      t_span4 = (3*tempo, 4*tempo)

      #comando para chamar o metodo de integracao
      sol4 = solve_ivp(ode_func4, t_span4, CI4, method=metodo, max_step = step)    # RK45 RK23 DOP853 BDF

      # EXTRAINDO OS RESULTADOS
      t4 = sol4.t
      m4 = sol4.y[0]
      T4 = sol4.y[1]
      S4 = sol4.y[2]
      X4 = sol4.y[3]

      #VISUALIZACO GRAFICA DOS 4 TEMPOS

      fig, axs = plt.subplots(2, 2, figsize=(12, 8))

      tempo_ms = tempo * 1000
      t1, t2, t3, t4 = t1 * 1000, t2 * 1000, t3 * 1000, t4 * 1000

      axs[0, 0].plot(t1, m1*1000, color='blue')
      axs[0, 0].plot(t2, m2*1000, color='blue')
      axs[0, 0].plot(t3, m3*1000, color='blue')
      axs[0, 0].plot(t4, m4*1000, color='blue')
      axs[0, 0].set_ylabel('Massa (g)')
      axs[0, 0].grid(True)

      axs[0, 1].plot(t1, T1-273.15, color='red')
      axs[0, 1].plot(t2, T2-273.15, color='red')
      axs[0, 1].plot(t3, T3-273.15, color='red')
      axs[0, 1].plot(t4, T4-273.15, color='red')
      axs[0, 1].set_ylabel('Temperatura (°C)')
      axs[0, 1].grid(True)

      axs[1, 0].plot(t1, S1, color='green')
      axs[1, 0].plot(t2, S2, color='green')
      axs[1, 0].plot(t3, S3, color='green')
      axs[1, 0].plot(t4, S4, color='green')
      axs[1, 0].set_ylabel('Entropia (J/kg·K)')
      axs[1, 0].grid(True)

      axs[1, 1].plot(t1, X1/1000, color='purple')
      axs[1, 1].plot(t2, X2/1000, color='purple')
      axs[1, 1].plot(t3, X3/1000, color='purple')
      axs[1, 1].plot(t4, X4/1000, color='purple')
      axs[1, 1].set_ylabel('Exergia (kJ)')
      axs[1, 1].grid(True)

      plt.tight_layout()

      buf = io.BytesIO()
      plt.savefig(buf, format="png")
      buf.seek(0)
      image_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
      buf.close()

      return render_template('simulador.html', plot=image_base64)
    else:
        return render_template('simulador.html')

if __name__ == '__main__':
    import os
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port)



