% the following matlab functions are used to calculate the Rand, 
% adjusted Rand, Wallace and other partition comparison coefficients
%
%    Copyright (C) 2009  UMMI@IMM
%
%    This file is part of Comparing Partitions website <http://www.comparingpartitions.info/>.    
%    Comparing Partitions website is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%[a,b,c,d,bc,dn,confmat,res]=PartAgreeCoef(c1,c2)
%Outputs:
%ri=rand Index
%AR = adjusted rand
%jac=jaccard
%w1- Wallace c1->c2
%w2 - Wallace c2->c1
%lar1 - Larsen c1->c2
%lar2 - Larsen c2->c1
%MH - Melia & Heckerman
%VI - variation of information
%nVI - normalized variation of information 
% Francisco Pinto fpinto@fm.ul.pt

function res=PartAgreeCoef_ARonly(c1,c2)

c1=c1-min(c1)+1;
c2=c2-min(c2)+1;
n=length(c1); % number cass  
ng1=max(c1);
ng2=max(c2);
%dn=n*(n-1)/2; %number of possible pairwise comparisons of cases(a+b+c+d)

%y1=pdist(c1(:),'ham');
%y2=pdist(c2(:),'ham');
%size(y1)
%size(y2)
%ad=sum((y1'==y2')); % number of pairwise concordances (matches (a) and mismatches(d))(a+d)
%bc=sum((y1'~=y2')); % number of pairwise discordances(b+c)
%Rand Index
%res.ri=ad/dn;

%a=sum((y1'==0).*(y2'==0));
%b=sum(y1'==0)-a;
%c=sum(y2'==0)-a;

%res.w1=a/sum(y1'==0);% check is =dn
%res.w1a=a/(a+b);
%res.w2a=a/(a+c);
%res.w2=a/sum(y2'==0);
%d=ad-a;

%Jaccard Index
%res.jac=a/(dn-d);

%confmat=crosstab(c1,c2);
confmat=full(sparse(c1,c2,1,ng1,ng2));

coltot=sum(confmat);
rowtot=sum(confmat')';
%summat=repmat(coltot,ng1,1)+repmat(rowtot,1,ng2);
%larsenmat=2*confmat./summat;
%res.lar1=mean(max((larsenmat'))');
%res.lar2=mean(max(larsenmat)');

%todelmat=larsenmat;
%cumval=0;
%for i=1:min([ng1;ng2])
%    [val]=max(max(todelmat)');
%    [rr,cc]=find(todelmat==val);
%        todelmat(rr(1),:)=0;
%        todelmat(:,cc(1))=0;
%        cumval=cumval+confmat(rr(1),cc(1));
%end
%res.MH=cumval/n;

%H1=-sum((rowtot/n).*log2((rowtot/n)));
%H2=-sum((coltot/n)'.*log2((coltot/n)'));
%indmat=(rowtot/n)*(coltot/n);
%nozeromat=(confmat/n)+(confmat==0);
%H12=-sum(sum((confmat/n).*log2(nozeromat)));
%MI=H1+H2-H12;
%res.VI=H1+H2-2*MI;
%res.nVI=res.VI/log2(n);

nis=sum(rowtot.^2);		%sum of squares of sums of rows
njs=sum(coltot.^2);		%sum of squares of sums of columns

t1=nchoosek(n,2);		%total number of pairs of entities
t2=sum(sum(confmat.^2));	%sum over rows & columnns of nij^2
t3=.5*(nis+njs);

%Expected index (for adjustment)
nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));

A=t1+t2-t3;		%no. agreements
%D=  -t2+t3;		%no. disagreements

if t1==nc
   res=0;			%avoid division by zero; if k=1, define Rand = 0
else
   res=(A-nc)/(t1-nc);		%adjusted Rand - Hubert & Arabie 1985
   %res.AR2=(ad-nc)/(dn-nc); % (a+d-nc)/(a+b+c+d-nc)
end
