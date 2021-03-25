
function es = esfuerzos(nodo1,nodo2,elemento)
    
    pos            = [nodo1.eqx ; nodo1.eqy ; nodo2.eqx ; nodo2.eqy];
    LargoOriginal  = sqrt( ((nodo1.x-nodo1.eqx)-(nodo2.x-nodo2.eqx))^2 + ((nodo1.y-nodo1.eqy)-(nodo2.y-nodo2.eqy))^2 ); %largo del elemento cuando esta sin esfuerzo
    LargoEsforzado = sqrt( (nodo1.x-nodo2.x)^2 + (nodo1.y-nodo2.y)^2 ); %Largo del elemento cuando soporta fuerzas
    
    %Si las fuerzas dan positivas significa que es tension, en cambio si da
    %negativa es compresion
    
    fuerzas  = elemento.k * (LargoEsforzado - LargoOriginal);
    es       = fuerzas;
end