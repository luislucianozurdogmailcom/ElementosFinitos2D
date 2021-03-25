%matrizRotacion
function matrizRotacion = Rotacion2D(phi)
matrizRotacion= [
    cos(phi)^2,cos(phi)*sin(phi),-cos(phi)^2,-cos(phi)*sin(phi);
    cos(phi)*sin(phi),sin(phi)^2,-cos(phi)*sin(phi),-sin(phi)^2;
    -cos(phi)^2,-cos(phi)*sin(phi),cos(phi)^2,cos(phi)*sin(phi);
    -cos(phi)*sin(phi),-sin(phi)^2,cos(phi)*sin(phi),sin(phi)^2];

s = size(matrizRotacion);
for i=1:s(1)
    for j=1:s(2)
        if(abs(matrizRotacion(i,j))<0.000000001)
            matrizRotacion(i,j) = 0;
        end
    end
end


end

