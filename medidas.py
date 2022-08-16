def delta(x,y):
    if(x==y):
        saida=1
    else:
        saida=0
    return saida

def auto(L,N):
    baseum=np.identity(N) # base canonica
    autovalores,autovetores=np.linalg.eig(L)
    #ordem crescente autovalores e autovetores
    idx=np.argsort(autovalores)
    autovalores=np.round(autovalores[idx],5)
    autovetores=np.round(autovetores[:,idx],5)
    baseum=baseum[:,idx]
    return baseum, autovetores,autovetores

# ortonormalizacao
for i in range(N):
	autovetores[:,i]=autovetores[:,i]/np.sqrt(np.sum(autovetores[:,i]**2))


    
def ineficienciamedia(autovalores,autovetores,base):
    INEFICIENCIA=0
    N1=len(autovalores)
    auxA=np.zeros(N)
    for a in range(N):
      aux=0
      for n in range(N):
        for m in range(N):
          aux=aux+(delta(autovalores[n],autovalores[m])*((sum(base[:,a]*autovetores[:,n])**2)*(sum(base[:,a]*autovetores[:,m])**2)))
      auxA[a]=aux
    INEFICIENCIA=sum(auxA)
    return INEFICIENCIA/N1

def Probabilidades_Medias (autovalores,autovetores,t,base,N):
  pmedioclassico=np.zeros(len(t))
  pmedioquanticoexato=np.zeros(len(t))
  pmedioquanticoaprox=np.zeros(len(t))
  aux1=np.zeros((N,N))
  for i in range(N):
    for k in range(N):
      aux1[i,k]=(sum(base[:,i]*autovetores[:,k]))**2
  
  for j in range(len(t)):

    pmedioclassico[j]=(1/N)*sum(np.exp(-(autovalores*t[j])))
    aux=0
    for a in range(N):
        aux=aux+((sum(np.cos(autovalores*t[j])*aux1[a,:])**2)+(sum(np.sin(autovalores*t[j])*aux1[a,:])**2))
    pmedioquanticoexato[j]=(1/N)*aux
    pmedioquanticoaprox[j]=(1/(N**2))*(((sum(np.sin(autovalores*t[j])))**2)+((sum(np.cos(autovalores*t[j])))**2))
  data=pd.DataFrame()
  data['t']=t
  data['pmedioclassico']=pmedioclassico
  data['pmedioquanticoaprox']=pmedioquanticoaprox
  data['pmedioquanticoexato']=pmedioquanticoexato
  return data

def Probabilidade_classica(autovalores,autovetores,base,N):
  t=10
  P=np.zeros((N,N))
  auxA=np.zeros((N,N))
  for a in range(N):
    for n in range(N):
      auxA[a,n]=(sum(base[:,a]*autovetores[:,n]))
  auxB=auxA
  for a in range(N):
    for b in range(N):
      P[a,b]=sum(np.exp(-(autovalores*t))*auxA[a,:]*auxB[b,:])
  return P

def Probabilidade_quantica(autovalores,autovetores,base,N):
  t=10
  PI=np.zeros((N,N))
  auxA=np.zeros((N,N))
  for a in range(N):
    for n in range(N):
      auxA[a,n]=(sum(base[:,a]*autovetores[:,n]))
  for a in range(N):
    for b in range(N):
      #PI[a,b]=((sum(np.cos(autovalores*t)*auxA[a,:]*auxB[b,:])**2)+(sum(np.sin(autovalores*t)*auxA[a,:]*auxB[b,:])**2))
       PI[a,b]=abs(sum(np.exp(autovalores*(-1j)*t)*auxA[a,:]*auxA[b,:]))**2 # Forma exponencial complexa

  return PI

def ineficiencia(autovalores,autovetores,base,N):
  INEF=np.zeros((N,N))

  auxA=np.zeros((N,N))
  for a in range(N):
    for n in range(N):
      auxA[a,n]=(sum(base[:,a]*autovetores[:,n]))

  for a in range(N):
    for b in range(N):
      aux=0
      for n in range(N):
        for m in range(N):
          aux=aux+(delta(autovalores[n],autovalores[m])*auxA[a,n]*auxA[b,n]*auxA[a,m]*auxA[b,m])
      INEF[a,b]=aux

  return INEF






    