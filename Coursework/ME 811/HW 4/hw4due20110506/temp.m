Ax = zeros(35,35);

for i = 1:5
	for j = 1:5
		J = j+5*(i-1)+5;
		I = J;
		Ax(J,I) = AxO(i,j);
		Ax(J,I-1) = AxW(i,j);
		Ax(J,I-5) = AxS(i,j);
		Ax(J,I+1) = AxE(i,j);
		Ax(J,I+5) = AxN(i,j);
	end
end
Ax = Ax(6:30,6:30);