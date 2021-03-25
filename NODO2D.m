classdef NODO2D
    properties
        x;
        y;
        eqx; %Posiciones respecto del equilibrio
        eqy;
    end
    methods
        function obj = NODO2D(a,b)
            %en relaciones se pone (posicion en x,posicion en y)
            obj.x   = a;
            obj.y   = b;
            obj.eqx = 0;
            obj.eqy = 0;
        end
    end
end