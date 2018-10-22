%loadsu2mesh

global ppp ttt  pointboun  filenmsu2  neumannBndryStr pointNeumboun
global mapNormalNeumBndry
pointboun=[];
pointNeumboun=[];

strele   ='NELEM=';
strpoin  ='NPOIN=';
strNmark ='NMARK=';
strMARKER_TAG  ='MARKER_TAG=';
strMARKER_ELEMS='MARKER_ELEMS=';
ne_line=0;
np_line=0;


nline=0;
bool1=0;
bool2=0;
bool3=0;

fid = fopen(filenmsu2, 'r');

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
        
        if strcmp(tline1(1:6),strNmark)
            n_Nmrksu2=str2num(strtrim(tline1(7:end)));
        %    np_line=nline;
            bool3=1;
        end        
        
    end
    
    if bool1==1 && bool2==1 && bool3==1
        break;
    end
end

fclose(fid);

cellMarkerTag=cell(n_Nmrksu2,1);
fid = fopen(filenmsu2, 'r');
jkmrk=1;
while feof(fid)==0 
    tline=fgetl(fid);
  %  nline=nline+1;
    tline1=strtrim(tline);
    if length(tline1)>11
        if strcmp(tline1(1:11),strMARKER_TAG)
            cellMarkerTag{jkmrk}=strtrim(tline1(12:end));
           % ne_line=nline;
           % bool1=1;
            jkmrk=jkmrk+1;
        end
    end
    
    if jkmrk> n_Nmrksu2
        break;
    end
end
fclose(fid);        

%cellMarkerElems=cell(n_Nmrksu2,1);
numBndryElems=zeros(n_Nmrksu2,1);

mapMarkerNumElems=containers.Map(cellMarkerTag(1:n_Nmrksu2),numBndryElems);
% for iii=1:n_Nmrksu2
%     mapMarkerNumElems(cellMarkerTag{iii})=0;
% end
cellBndryElems=cell(n_Nmrksu2,1);
% for i=1:n_Nmrksu2
%     cellBndryElems{i}=[];
% end
mapBndryElems=containers.Map(cellMarkerTag(1:n_Nmrksu2),cellBndryElems); 
% every kind of boundary, the boundary edge, here Elems are same as the su2
% for boundary edge in 2d
% cell to save the key string, should use mapBndryElems(cellMarkerTag{1})to
% get the value(not mapBndryElems(cellMarkerTag(1))  )

%mapBndryNormal=containers.Map(cellMarkerTag,cell(n_Nmrksu2,1));

fid = fopen(filenmsu2, 'r');
jkmrk=1;
while feof(fid)==0
    tline=fgetl(fid);
    %  nline=nline+1;
    tline1=strtrim(tline);
    if length(tline1)>11
        if strcmp(tline1(1:11),strMARKER_TAG)
            strTmpLine=strtrim(tline1(12:end));
            bltmp1=1;
            
            while bltmp1==1
                tline=fgetl(fid);
                %  nline=nline+1;
                tline1=strtrim(tline);
                if strcmp(tline1(1:13),strMARKER_ELEMS)
                    mapMarkerNumElems(strTmpLine)=str2num(strtrim(tline1(14:end)));
                    bltmp1=0;
                end
            end
            
            jkmrk=jkmrk+1;
        end
    end
             
    if jkmrk> n_Nmrksu2
        break;
    end
end
fclose(fid);  

fid = fopen(filenmsu2, 'r');
jkmrk=1;
while feof(fid)==0
    tline=fgetl(fid);
    %  nline=nline+1;
    tline1=strtrim(tline);
    if length(tline1)>11
        if strcmp(tline1(1:11),strMARKER_TAG)
            strTmpLine=strtrim(tline1(12:end));
            jkmrk=jkmrk+1;
            tline=fgetl(fid);
            bltmp1=1;
            arr4BndryEdges=zeros(mapMarkerNumElems(strTmpLine),2); % for 2d
            
            while bltmp1<=mapMarkerNumElems(strTmpLine)
                tline=fgetl(fid);
                %  nline=nline+1;
                tline1=strtrim(tline);
                arrTmp=str2num(tline1);
                arr4BndryEdges(bltmp1,:)=  arrTmp(2:3);             
                bltmp1=bltmp1+1;
            end  
            
            mapBndryElems(strTmpLine)=arr4BndryEdges;
        end
    end
             
    if jkmrk> n_Nmrksu2
        break;
    end
end
fclose(fid);  

%%%%%%%%%%%%%%%%%%%%%----------
pts_e=zeros(n_elsu2,5);
pts=zeros(n_posu2,3);
nline=0;
ij=1;
ij2=1;
fid = fopen(filenmsu2, 'r');
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
cd rbfdq;
neb=size(eboun,1);
for i=1:neb
    pointboun=[pointboun eboun(i,:)];
end
pointboun=sort(pointboun);
pointboun=pointboun';
pointboun=unique(pointboun);


if isKey(mapBndryElems,neumannBndryStr)
    aaatm=mapBndryElems(neumannBndryStr)+1;
    pointNeumboun=sort(unique(aaatm(:)));
    mapNormalNeumBndry=containers.Map(pointNeumboun,cell(length(pointNeumboun),1));
    for ils1=1:length(pointNeumboun)
        mapNormalNeumBndry(pointNeumboun(ils1))=zeros(1,2);
    end
    
    mapNumEdgeBP=containers.Map(pointNeumboun,zeros(length(pointNeumboun),1)); 
    for itmp=1:size(aaatm,1)
        dy=ppp(aaatm(itmp,2),2)-ppp(aaatm(itmp,1),2);
        dx=ppp(aaatm(itmp,2),1)-ppp(aaatm(itmp,1),1);
        r2dxdy=sqrt(dx*dx+dy*dy);
        mapNumEdgeBP(aaatm(itmp,2))=mapNumEdgeBP(aaatm(itmp,2))+1;
        mapNumEdgeBP(aaatm(itmp,1))=mapNumEdgeBP(aaatm(itmp,1))+1;
        
%         mapNormalNeumBndry(aaatm(itmp,2))= ...
%             (mapNormalNeumBndry(aaatm(itmp,2))+[dy, -dx]/r2dxdy)/mapNumEdgeBP(aaatm(itmp,2));
%         mapNormalNeumBndry(aaatm(itmp,1))= ...
%             (mapNormalNeumBndry(aaatm(itmp,1))+[dy, -dx]/r2dxdy)/mapNumEdgeBP(aaatm(itmp,1));
       
       if mapNumEdgeBP(aaatm(itmp,2))==1
           mapNormalNeumBndry(aaatm(itmp,2))=[dy, -dx]/r2dxdy;
       end
       if mapNumEdgeBP(aaatm(itmp,1))==1
           mapNormalNeumBndry(aaatm(itmp,1))=[dy, -dx]/r2dxdy;
       end
      if  mapNumEdgeBP(aaatm(itmp,2))==2
        mapNormalNeumBndry(aaatm(itmp,2))= ...
            (mapNormalNeumBndry(aaatm(itmp,2))+[dy, -dx]/r2dxdy)/mapNumEdgeBP(aaatm(itmp,2));
      end
      if  mapNumEdgeBP(aaatm(itmp,1))==2      
      mapNormalNeumBndry(aaatm(itmp,1))= ...
            (mapNormalNeumBndry(aaatm(itmp,1))+[dy, -dx]/r2dxdy)/mapNumEdgeBP(aaatm(itmp,1));    
      end
        
    end
    
end
 
for itmp=1:length(pointNeumboun)
    r2dxdy=mapNormalNeumBndry(pointNeumboun(itmp));
    mapNormalNeumBndry(pointNeumboun(itmp))=r2dxdy/norm(r2dxdy);
end
    


clear eboun;
fprintf('------(Done mesh generation.)-------\n\n')


    