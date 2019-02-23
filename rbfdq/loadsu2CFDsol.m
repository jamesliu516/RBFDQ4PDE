%loadsu2CFDsol.m su2 v6.1

global ppp  filenmsu2Sol vecVel
np1=size(ppp,1);
vecVel=zeros(np1,2);

fid = fopen(filenmsu2Sol, 'r');
nline=0;

while feof(fid)==0 
    tline=fgetl(fid);
    nline=nline+1;
    tline1=strtrim(tline);
    
    if nline >= 4 && nline<= np1+3
       linedata=str2num(tline1);
       vecVel(nline-3,:)=linedata(4:5);
    end
    
    if nline >np1+3
        break;
    end
end
fclose(fid);



fprintf('+         SU2v6.1 Done CFD solution loading.                +\n');


    