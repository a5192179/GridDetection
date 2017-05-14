tic
for i=1:10
    a=load('a.txt');
end
toc
tic
for i=1:10
    a=importdata('a.txt');
end
toc