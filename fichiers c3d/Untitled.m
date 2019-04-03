[X,Y,Z] = sphere(16);
x = [0.5*X(:); 0.75*X(:); X(:)];
y = [0.5*Y(:); 0.75*Y(:); Y(:)];
z = [0.5*Z(:); 0.75*Z(:); Z(:)];



S = repmat([100,50,5],numel(X),1);
s = S(:);



figure
scatter3(x,y,z,s)
view(40,35)