# lab4

REPORTE LABORATORIO 04

INTRODUCCIÓN 	

Los tests estadísticos se utilizan para validar o rechazar las hipótesis de modelación. Tratan de distinguir entre lo qué brinda un resultado satisfactorio de lo qué brinda un resultado erróneo en el marco del modelo dado.  Existen distintos tipos de tests estadísticos, en donde se determinan diferentes requisitos de los modelos creados. 

En este laboratorio se implementaron diez tests con el objetivo de evaluar la calidad de los generadores pseudoaleatorios creados en el laboratorio pasado. Luego se muestra una tabla con el resumen de los tests estadísticos y por último una gráfica en donde se pueden observar los tests qué fallaron. Esto se realizó con los tres  ejemplos de cadenas qué se implementaron en el Laboratorio 3. 

METODOLOGÍA

Batería de diez tests estadísticos 
Se implementaron diez tests, los cuales se encuentran en los links de github brindado por el licenciado, a cada uno se le realizaron las modificaciones necesarias para funcionar con los generadores pseudoaleatorios. Los tests implementados son: 
monobit_test.py
frequency_within_block_test.py
runs_test.py
longest_run_ones_in_a_block_test.py
binary_matrix_rank_test.py
linear_complexity_test.py
random_excursion_variant_test.py
random_excursion_test.py
approximate_entropy_test.py
random_excursion_test.py

Tabla de resumen de los diez tests estadísticos
Se realizó una función que, dada una cadena de bits, elabora una tabla con el resumen de los diez tests estadísticos implementados, esta función se aplicó a los ejemplos del Laboratorio 03 para poder tener una métrica de evaluación de los mismos. 

Histograma de frecuencias 
Se realizó una función qué con los generadores pseudo-aleatorios implementados en el Laboratorio 03, genera mil cadenas pseudo-aleatorias diversas y se ejecuta la batería de test para estas cadenas, luego se grafica un histograma de frecuencias del número de tests fallados por las cadenas aleatorias.


RESULTADOS 

LCG

Corrida con 1000 iteraciones


Se obtuvo la tabla de resumen aplicada a una cadena aleatoria utilizando el generador LCG y como se puede observar paso los 10 tests. En cuanto al histograma se observa que solo la prueba 8 no dio ningún error y todas las demás si, por lo que tomando en cuenta esto podríamos decir que no es muy confiable este generador. Sin embargo, observando en el histograma las frecuencias de error de cada prueba podemos notar que son bajas, llegando a ser las más altas la prueba 7 con 67 y la prueba 10 con 83 respectivamente, por lo que podría llegar a considerarse como un generador confiable.

Wichman Hill

Corrida con 1000 iteraciones






Corrida con 5 en vez de 1000 iteraciones

LFSR
Tablas de resumen 

	Cadena Pequeña 
	
	
Cadena Mediana 
	

Cadena Larga 
	

Histograma de Frecuencias 
	
	
Con las mil cadenas aleatorias se obtuvieron resultados qué nos indica qué el generador pseudoaleatorio no es muy eficiente con todos los tests ya qué de los diez qué se implementaron seis dieron error siendo el décimo con una mayor frecuencia. Sin embargo se realizaron pruebas de los test en donde se utilizaron cadenas de utilizaron 2000000 bits y todos dieron resultados positivos, por lo cual no se tiene un histograma de frecuencias pero si la tabla de resumen. 


CONCLUSIONES
   
  
El generador pseudoaleatorio LFSR implementado es más eficiente cuando se trabaja con cadenas muy largas de bits y al igual qué los resultados obtenidos cifrando las imágenes el en Laboratorio 03, podemos observar qué mientras más larga sea la cadena, mejor son los resultados de los tests. 

LINK DEL REPOSITORIO
https://github.com/JosePabloPonce/lab4

