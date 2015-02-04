test.Ax.O = abs(Chris.links.mom.o-Ax.O) > 100*eps;
test.Ax.E = abs(Chris.links.mom.e-Ax.E) > 100*eps;
test.Ax.W = abs(Chris.links.mom.w-Ax.W) > 100*eps;
test.Ax.N = abs(Chris.links.mom.n-Ax.N) > 100*eps;
test.Ax.S = abs(Chris.links.mom.s-Ax.S) > 100*eps;
test.Ax.P.x = abs(Chris.links.mom.p.x-Ax.P) > 100*eps;
test.Ax.P.y = abs(Chris.links.mom.p.y-Ay.P) > 100*eps;

test.Ap.O = abs(Chris.links.p.o-Ap.O) > 100*eps;
test.Ap.E = abs(Chris.links.p.e-Ap.E) > 100*eps;
test.Ap.W = abs(Chris.links.p.w-Ap.W) > 100*eps;
test.Ap.N = abs(Chris.links.p.n-Ap.N) > 100*eps;
test.Ap.S = abs(Chris.links.p.s-Ap.S) > 100*eps;
test.Ap.P = abs(Chris.links.p.p-Ap.P) > 100*eps;

test.At.O = abs(Chris.links.T.o-At.O) > 100*eps;
test.At.E = abs(Chris.links.T.e-At.E) > 100*eps;
test.At.W = abs(Chris.links.T.w-At.W) > 100*eps;
test.At.N = abs(Chris.links.T.n-At.N) > 100*eps;
test.At.S = abs(Chris.links.T.s-At.S) > 100*eps;
test.At.P = abs(Chris.links.T.p-At.P) > 100*eps;

test.u = abs(Chris.u.c-u) > 100*eps;
test.v = abs(Chris.v.c-v) > 100*eps;
test.t = abs(Chris.T.c-T) > 100*eps;

test.uf.e = abs(Chris.u.f(:,2:end)-uf.e) > 100*eps;
test.uf.w = abs(Chris.u.f(:,1:end-1)-uf.w) > 100*eps;
test.vf.n = abs(Chris.v.f(2:end,:) -vf.n) > 100*eps;
test.vf.s = abs(Chris.v.f(1:end-1,:) -vf.s) > 100*eps;
test.p = abs(Chris.p.c-p) > 100*eps;


summation.Ax.O = sum(sum(test.Ax.O));
summation.Ax.W = sum(sum(test.Ax.W));
summation.Ax.E = sum(sum(test.Ax.E));
summation.Ax.N = sum(sum(test.Ax.N));
summation.Ax.S = sum(sum(test.Ax.S));
summation.Ax.P.x = sum(sum(test.Ax.P.x));
summation.Ax.P.y = sum(sum(test.Ax.P.y));

summation.Ap.O = sum(sum(test.Ap.O));
summation.Ap.W = sum(sum(test.Ap.W));
summation.Ap.E = sum(sum(test.Ap.E));
summation.Ap.N = sum(sum(test.Ap.N));
summation.Ap.S = sum(sum(test.Ap.S));
summation.Ap.P = sum(sum(test.Ap.P));

summation.At.O = sum(sum(test.At.O));
summation.At.W = sum(sum(test.At.W));
summation.At.E = sum(sum(test.At.E));
summation.At.N = sum(sum(test.At.N));
summation.At.S = sum(sum(test.At.S));
summation.At.P = sum(sum(test.At.P));

summation.u = sum(sum(test.u));
summation.v = sum(sum(test.v));
summation.t = sum(sum(test.t));

summation.uf.e = sum(sum(test.uf.e));
summation.uf.w = sum(sum(test.uf.w));
summation.vf.n = sum(sum(test.vf.n));
summation.vf.s = sum(sum(test.vf.s));
summation.p = sum(sum(test.p));
