# lab4	

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

