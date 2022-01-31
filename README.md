# ExamplesFCCee

Ejemplos de código para empezar a trabajar con sucesos e+e-\rightarrow ZH  (FCCee) 

El objetivo es seguir este estudio sobre el recoil del Higgs: https://github.com/HEP-FCC/FCCeePhysicsPerformance/tree/ZH_recoil/case-studies/higgs/mH-recoil 

Vamos a reproducir los plots del recoil y de la masa del Z a partir de las muestras FCCee 'reducidas' creadas por Juan Alcaraz y almacenadas en las maquinas del CIEMAT.

###  Como correr los ejemplos

Desde una maquina ciemat conectate a las user interfaces 
```
ssh gaeui01   
```

Descarga el codigo de este tutorial

```
git clone https://github.com/mcepeda/ExamplesFCCee.git 
cd ExamplesFCCee/
```

Carga root
```
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-centos7-gcc48-opt/setup.sh
```

Corre los ejemplos:

a) Ejemplo con loops: 
```
python -i FirstSteps/plot_example_loop.py   
```

b) Ejemplo con dataframes:
```
python -i FirstSteps/plot_example_dataframe.py
```


