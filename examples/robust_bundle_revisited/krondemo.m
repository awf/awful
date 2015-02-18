%function krondemo


m=21;
n=73;
r = 3;

profile on

A = rand(m*r,n);
B = rand(n, n*r);

H = kron(A,B);

whos A B H

krondemo_test(A,B,H)

profile viewer
