%Programa para simular una viga larga

Longitud            = input('Introduzca la longitud de la viga : ');
Altitud             = input('Introduzca la altitud de la viga : ');  
BloquesVerticales   = input('Introduzca la cantidad de celdas verticales: ');
BloquesHorizontales = input('Introduzca la cantidad de celdas horizontales: ');

%Deje predefinido un modulo de Young y un area transverzal

AltitudCelda    = Altitud / BloquesVerticales;
LongitudCelda   = Longitud / BloquesHorizontales;
AreaTransversal = 0.01;
E               = 1.5e11;
kVertical       = (AreaTransversal*E)/AltitudCelda;
kHorizontal     = (AreaTransversal*E)/LongitudCelda;
kDiagonal       = (AreaTransversal*E)/(sqrt( AltitudCelda^2 + LongitudCelda^2 ));
anguloDiagonal  = atan(AltitudCelda / LongitudCelda);

%Creo la matriz Global y calculo constantes del sistema, como nodos,elementos, etc

Bloques                       = BloquesHorizontales * BloquesVerticales;
CantidadNodos                 = ((BloquesVerticales + 1)*(BloquesHorizontales + 1));
MatrizGlobal                  = MatrizGlobal( CantidadNodos , 2);
CantidadElementosVerticales   = BloquesVerticales * (BloquesHorizontales + 1);
CantidadElementosHorizontales = BloquesHorizontales * (BloquesVerticales + 1);
CantidadElementosDiagonales   = Bloques * 2;

%Creo todos los Elementos verticales:

for i=1:CantidadElementosVerticales
    if i==1
        El = ELEMENTO(1,((BloquesVerticales - 1)*(BloquesHorizontales + 1) + 3),kVertical,pi / 2);
    elseif i==2 
        El = [El ELEMENTO(2,((BloquesVerticales - 1)*(BloquesHorizontales + 1) + 3 + BloquesHorizontales),kVertical,pi / 2)];
    elseif i <= CantidadElementosVerticales - (BloquesHorizontales + 1 - 2)
        El = [El ELEMENTO(i,i+BloquesHorizontales+1,kVertical,pi/2)]; 
    else
        El = [El ELEMENTO(i+1, i+BloquesHorizontales+1, kVertical,pi/2)];
    end
end

%Creo todos los Elementos horizontales:

for i=1:CantidadNodos
    if i==1
        El = [El ELEMENTO(1,((BloquesVerticales)*(BloquesHorizontales + 1) + 3),kHorizontal,0)];
    elseif i==2 
        El = [El ELEMENTO(2,CantidadNodos,kHorizontal,0)];
    elseif i < CantidadNodos - BloquesHorizontales + 1
        
        if mod(i-2,BloquesHorizontales+1) == 0 
            %No hace nada
        else
            El = [El ELEMENTO(i ,i+1 ,kHorizontal,0)];
        end
    
    elseif i >= CantidadNodos - BloquesHorizontales + 2 && i ~= CantidadNodos
        El = [El ELEMENTO(i, i+1, kHorizontal,0)];
    end
end

%Creo todos los Elementos a "45�":

for i=1:CantidadNodos
    if i==1
        El = [El ELEMENTO(1,((BloquesVerticales - 1)*(BloquesHorizontales + 1) + 4),kDiagonal,anguloDiagonal)];
    elseif i >= BloquesHorizontales + 1 + 3 && i < (BloquesHorizontales + 1)*BloquesVerticales + 3

        if mod(i-2,BloquesHorizontales+1) == 0 
            %No hace nada
        else
            El = [El ELEMENTO(i ,i - BloquesHorizontales ,kDiagonal,anguloDiagonal)];
        end
    
    elseif (i >= (BloquesHorizontales + 1) * ( BloquesVerticales) + 3)
        El = [El ELEMENTO(i, i - BloquesHorizontales + 1, kDiagonal,anguloDiagonal)];
    end
end

%Creo todos los Elementos a "-45�":

for i=1:CantidadNodos
    if i==2
        El = [El ELEMENTO(2,CantidadNodos - BloquesHorizontales ,kDiagonal,-anguloDiagonal)];
    elseif i >= 3 && i < (BloquesHorizontales+1)*(BloquesVerticales - 1) + 3

        if mod(i-2,BloquesHorizontales+1) == 0 
            %No hace nada
        else
            El = [El ELEMENTO(i ,i + BloquesHorizontales + 2 ,kDiagonal, -anguloDiagonal)];
        end
    
    elseif i >= (BloquesHorizontales+1)*(BloquesVerticales - 1) + 3 && i < (BloquesHorizontales+1)*BloquesVerticales + 1
        El = [El ELEMENTO(i, i + BloquesHorizontales + 1, kDiagonal, -anguloDiagonal)];
    end 
end

%Creo Los nodos para cada punto de la viga que se discretiz�

%para el grafico de nodos
varx  = [];
vary  = [];
Nodos = [];
Const = 0;
for i=1:CantidadNodos
    if(i==1)
        Nodos = [Nodos NODO2D(0,0)];
    
    elseif (i==2)
        Nodos = [Nodos NODO2D(Longitud,0)];
    
    elseif (i>=3 && i <= (BloquesHorizontales + 1)*(BloquesVerticales)+2 )
        if(mod(i-3,BloquesHorizontales+1) == 0 && i ~= 3)
            Const = 1 + Const;
            Nodos = [Nodos NODO2D((i - 3 - Const*(BloquesHorizontales+1))*LongitudCelda, Altitud - Const*AltitudCelda)];
        else
            Nodos = [Nodos NODO2D((i - 3 - Const*(BloquesHorizontales+1))*LongitudCelda, Altitud - Const*AltitudCelda)];
        end
    else
        Nodos = [Nodos NODO2D((i - 2 - (Const+1)*(BloquesHorizontales+1))*LongitudCelda,0)];
    end
    varx = [varx Nodos(i).x]; %para el grafico de nodos
    vary = [vary Nodos(i).y]; %para el grafico de nodos
end



%Ensamblo la matriz global 

[s1,s2] = size(El);
for i=1:s2
    MatrizGlobal=insercion(MatrizGlobal,El(i).matrizLocal,El(i).nodos(1),El(i).nodos(2));
end

%Armo la submatriz KF

s3 = size(MatrizGlobal);
KF = MatrizGlobal(5:s3(1),5:s3(1));

%Armo la submatriz KEF

KEF = MatrizGlobal(1:4,5:s3(1));

%Creo la matriz columna de fuerzas conocidas

F      = [];
fuerza = input('Por favor introduzca la fuerza soportada :');

% -------------Grafico para ubicar visualmente los nodos ------------
figure(1)
plot(varx,vary,"color","none","marker","o","markerFaceColor","k"); %Grafico para identificar el nodo
hold on
for i=1:CantidadNodos
    text(Nodos(i).x,Nodos(i).y,['\leftarrow' num2str(i)]);
end
%xlim([-LongitudViga/25 LongitudViga+LongitudViga/25])
%ylim([-AltitudViga/5 AltitudViga+AltitudViga/5])
hold off
%------------------------------------------------------------------

choice = menu('Elija una forma de carga','Puntual','Distribuida');

if (choice == 1)
    Nodo   = input('Por favor introduzca el nodo que la soporta :');
    dim    = input('Por favor introduzca la dimension en la que esta la fuerza siendo 1:x y 2:y :');

    for i=1:CantidadNodos*2-4
        if (i==Nodo*2-6+dim) 
            F = [F;fuerza];
        else
            F = [F;0];
        end
    end
    
else
    carga = fuerza/CantidadNodos;
    
    for i=1:CantidadNodos*2-4
        if (mod(i,2) == 0 ) 
            F = [F;carga];
        else
            F = [F;0];
        end
    end
end

%Creo el sistema y lo resuelvo para posiciones desconocidas
sistema    = [KF F];
sistema    = rref(sistema);
tam        = size(sistema);
posiciones = sistema(:,tam(2));

%Luego utilizo las posiciones encontradas guardadas en posiciones para usarla
%para encontrar las reacciones con la submatriz KEF
reacciones = KEF * posiciones;

for i=2:2:4
    fprintf("Las reacciones para el nodo %d es:(%d,%d) \n",i/2,reacciones(i-1),reacciones(i));
end

for i=2:2:CantidadNodos*2-4
    fprintf("Las posiciones respecto del punto de equilibro para el nodo %d son:(%d,%d) \n",i/2+2,posiciones(i-1),posiciones(i));
    Nodos((i/2)+2).x   = Nodos((i/2)+2).x + posiciones(i-1);
    Nodos((i/2)+2).eqx = posiciones(i-1);
    Nodos((i/2)+2).y   = Nodos((i/2)+2).y + posiciones(i);
    Nodos((i/2)+2).eqy = posiciones(i);
    
end

%Calculo los esfuerzos para los elementos

for i=1:s2
    %Uso la funcion esfuerzos, poniendo en ella los nodos relacionados
    %y el elemento.
    El(i).esfuerzo = esfuerzos(Nodos(El(i).nodos(1)) , Nodos(El(i).nodos(2)) , El(i));
    
end

%Seccion de ploteo del sistema 

matCx = [0;0]; matCy = [0;0];
matTx = [0;0]; matTy = [0;0];
for i=1:s2
    if(El(i).esfuerzo < 0)
        a     = [Nodos(El(i).nodos(1)).x ; Nodos(El(i).nodos(2)).x];
        matCx = [matCx a];
        a     = [Nodos(El(i).nodos(1)).y ; Nodos(El(i).nodos(2)).y];
        matCy = [matCy a];
    else
        a     = [Nodos(El(i).nodos(1)).x ; Nodos(El(i).nodos(2)).x];
        matTx = [matTx , a];
        a     = [Nodos(El(i).nodos(1)).y ; Nodos(El(i).nodos(2)).y];
        matTy = [matTy , a];
    end
end

figure(2)
plot(matCx,matCy,'r',matTx,matTy,'b');
hold on
plot(matCx,matCy,'color','none','marker','o','markerFaceColor','k','markerSize',5);
plot(matTx,matTy,'color','none','marker','o','markerFaceColor','k','markerSize',5);


for i=1:CantidadNodos
    text(Nodos(i).x,Nodos(i).y,['\leftarrow' num2str(i)]);
end
xlim([-Longitud/25 Longitud+Longitud/25])
ylim([-Altitud/25 Altitud+Altitud/25])

hold off