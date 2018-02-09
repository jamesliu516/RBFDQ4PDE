%loadsu2mesh

global ppp ttt  pointboun
pointboun=[];
filenm='circleu1sss.su2';
fid = fopen(filenm, 'r');
strele='NELEM=';
strpoin='NPOIN=';
ne_line=0;
np_line=0;
nline=0;
bool1=0;
bool2=0;
while feof(fid)==0 
    tline=fgetl(fid);
    nline=nline+1;
    tline1=strtrim(tline);
    if length(tline1)>6
        if strcmp(tline1(1:6),strele)
            n_elsu2=str2num(strtrim(tline1(7:end)));
            ne_line=nline;
            bool1=1;
        end
        
        if strcmp(tline1(1:6),strpoin)
            n_posu2=str2num(strtrim(tline1(7:end)));
            np_line=nline;
            bool2=1;
        end
    end
    
    if bool1==1 && bool2==1
        break;
    end
end
fclose(fid);

pts_e=zeros(n_elsu2,5);
pts=zeros(n_posu2,3);
nline=0;
ij=1;
ij2=1;
fid = fopen(filenm, 'r');
while feof(fid)==0 
    tline=fgetl(fid);
    tline=strtrim(tline);
    nline=nline+1;
    if nline>=ne_line+1 && nline<=ne_line+n_elsu2
        pts_e(ij,:)=str2num(tline);
        ij=ij+1;
    end
    
    if nline>=np_line+1 && nline<=np_line+n_posu2
        pts(ij2,:)=str2num(tline);
        ij2=ij2+1;
    end    
end

fclose(fid);

% ppp=zeros(n_posu2,2);
% ttt=zeros(n_elsu2,3);
ppp=pts(:,1:2);
ttt=pts_e(:,2:4)+1;

cd ../
eboun=boundedges(ppp,ttt);
neb=size(eboun,1);
for i=1:neb
    pointboun=[pointboun eboun(i,:)];
end
pointboun=sort(pointboun);
pointboun=pointboun';
pointboun=unique(pointboun);

clear eboun;
fprintf('(Done mesh generation.)\n\n')
cd rbfdq

    