% Funcion insertar en la global;
function matrizG = insercion(matrizGlobal,A,nodo1,nodo2)

tamanos = size(matrizGlobal);

for i=1:tamanos(1)
    for j=1:tamanos(2)
        if(i == 2*nodo1-1 && j == 2*nodo1-1) % para 1,1 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(1,1);
        end
        
        if (i == 2*nodo1-1 && j == 2*nodo1)% para 1,2 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(1,2);
        end
        
        if (i == 2*nodo1 && j == 2*nodo1-1)% para 2,1 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(2,1);
        end
        
        if (i == 2*nodo1 && j == 2*nodo1)% para 2,2 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(2,2);
        end
        
        if(i == 2*nodo1-1 && j == 2*nodo2-1)% para 1,3 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(1,3);
        end
        
        if(i == 2*nodo1-1 && j == 2*nodo2)% para 1,4 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(1,4);
        end
        
        if(i == 2*nodo1 && j == 2*nodo2-1)% para 2,3 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(2,3);
        end
        
        if (i == 2*nodo1 && j == 2*nodo2)% para 2,4 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(2,4);
        end
        
        if (i == 2*nodo2-1 && j == 2*nodo1-1)% para 3,1 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(3,1);
        end
        
        if (i == 2*nodo2-1 && j == 2*nodo1)% para 3,2 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(3,2);
        end
        
        if (i == 2*nodo2 && j == 2*nodo1-1)% para 4,1 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(4,1);
        end
        
        if (i == 2*nodo2 && j == 2*nodo1)% para 4,2 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(4,2);
        end
        
        if (i == 2*nodo2-1 && j == 2*nodo2-1)% para 3,3 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(3,3);
        end
        
        if (i == 2*nodo2-1 && j == 2*nodo2)% para 3,4 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(3,4);
        end
        
        if (i == 2*nodo2 && j == 2*nodo2-1)% para 4,3 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(4,3);
        end
        
        if (i == 2*nodo2 && j == 2*nodo2)% para 4,4 de la matriz local
            matrizGlobal(i,j) = matrizGlobal(i,j) + A(4,4);
        end
    end
end

matrizG = matrizGlobal;

end