classdef ELEMENTO
    properties
        nodos
        k
        angulo
        matrizLocal
        esfuerzo
    end
    methods
        function obj = ELEMENTO(a,b,k,angulo)
            obj.nodos       = [a,b];
            obj.k           = k;
            obj.angulo      = angulo;
            obj.matrizLocal = obj.k*Rotacion2D(angulo);
            obj.esfuerzo    = 0;
        end
    end
end