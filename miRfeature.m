function feature=miRfeature(SS)
n=length(SS);
feature=zeros(n,42);
for i=1:n
  A=cell2mat(cellstr(SS(i,3)));
  p1=findstr(SS(i,3),'(');
  p2=findstr(SS(i,3),')');
  basepair=length(p1);
  Len=length(A);
  CenterLoop=p2(1)-p1(basepair)-1;
  arm3Len=Len-CenterLoop-p1(basepair);
  SD=abs(arm3Len-p1(basepair));
  LBR=Len/basepair;
  Seq=char(SS(i,2));
  base=basecount(Seq);
  [b,MFE]=rnafold(Seq);
  site=findstr(SS(i,3),'.');
  SL=length(site);
  if site(1)==1 && site(SL)==Len
      tailcount=2;
  elseif site(1)==1 && site(SL)~=Len
      tailcount=1;
  elseif site(1)~=1 && site(SL)==Len
      tailcount=1;
  elseif site(1)~=1 && site(SL)~=Len
      tailcount=0;
  end
  if tailcount==0
      taillength=0;
  else
      tp1=split(SS(i,3),'(');
      tp2=split(SS(i,3),')');
      Tp1=char(tp1(1));
      tailLen1=length(Tp1);
      Tp2=char(tp2(length(tp2)));
      tailLen2=length(Tp2); 
      if tailLen1>=tailLen2
          taillength=tailLen1;
      else
          taillength=tailLen2;
      end
  end
  AAA=char(SS(i,3));
  b1=AAA(1:p1(length(p1)));
  b2=AAA(p2(1):Len);
  bb1=split(b1,'.');
  bb2=split(b2,'.');
  bb1=string(bb1);
  bb2=string(bb2);
  bulgecount1=length(find(bb1~=''))-1;
  bulgecount2=length(find(bb2~=''))-1;
  if bulgecount1>=bulgecount2
      bulgercount=bulgecount1;
  else
      bulgercount=bulgecount2;
  end
  feature(i,1)=SD;
  feature(i,2)=basepair;
  feature(i,3)=(base.C+base.G)/(base.A+base.C+base.G+base.T);%%G+C
  feature(i,4)=LBR;
  feature(i,5)=Len;
  feature(i,6)=CenterLoop;
  feature(i,7)=MFE/Len;%%ÄÜ»ù±È
  feature(i,8)=tailcount;
  feature(i,9)=taillength;
  feature(i,10)= bulgercount;
end
for k=1:n
    b=char(SS(k,3));
    seq=char(SS(k,2));
    Len=length(seq);
    for j=1:Len
        if b(j)==')'
            b(j)='(';
        end
    end
    count=zeros(1,32);
    for i=1:Len-2
        if b(i:i+2)=='(((' & seq(i+1)=='A'
            count(1)=count(1)+1;
        elseif b(i:i+2)=='(((' & seq(i+1)=='U'
            count(2)=count(2)+1;
        elseif b(i:i+2)=='(((' & seq(i+1)=='C'
            count(3)=count(3)+1;
        elseif b(i:i+2)=='(((' & seq(i+1)=='G'
            count(4)=count(4)+1;
        elseif b(i:i+2)=='((.' & seq(i+1)=='A'
             count(5)=count(5)+1;
        elseif b(i:i+2)=='((.' & seq(i+1)=='U'
             count(6)=count(6)+1;
        elseif b(i:i+2)=='((.' & seq(i+1)=='C'
            count(7)=count(7)+1;
        elseif b(i:i+2)=='((.' & seq(i+1)=='G'
            count(8)=count(8)+1;
        elseif b(i:i+2)=='(..' & seq(i+1)=='A'
            count(9)=count(9)+1;
        elseif b(i:i+2)=='(..' & seq(i+1)=='U'
            count(10)=count(10)+1;
        elseif b(i:i+2)=='(..' & seq(i+1)=='C'
            count(11)=count(11)+1;
        elseif b(i:i+2)=='(..' & seq(i+1)=='G'
            count(12)=count(12)+1;
        elseif b(i:i+2)=='(.(' & seq(i+1)=='A'
            count(13)=count(13)+1;
        elseif b(i:i+2)=='(.(' & seq(i+1)=='U'
            count(14)=count(14)+1;
        elseif b(i:i+2)=='(.(' & seq(i+1)=='C'
            count(15)=count(15)+1;
        elseif b(i:i+2)=='(.(' & seq(i+1)=='G'
            count(16)=count(16)+1;
        elseif b(i:i+2)=='...' & seq(i+1)=='A'
            count(17)=count(17)+1;
        elseif b(i:i+2)=='...' & seq(i+1)=='U'
            count(18)=count(18)+1;
        elseif b(i:i+2)=='...' & seq(i+1)=='C'
            count(19)=count(19)+1;
        elseif b(i:i+2)=='...' & seq(i+1)=='G'
            count(20)=count(20)+1;
        elseif b(i:i+2)=='..(' & seq(i+1)=='A'
            count(21)=count(21)+1;
        elseif b(i:i+2)=='..(' & seq(i+1)=='U'
            count(22)=count(22)+1;
        elseif b(i:i+2)=='..(' & seq(i+1)=='C'
            count(23)=count(23)+1;
        elseif b(i:i+2)=='..(' & seq(i+1)=='G'
            count(24)=count(24)+1;
        elseif b(i:i+2)=='.(.' & seq(i+1)=='A'
            count(25)=count(25)+1;
        elseif b(i:i+2)=='.(.' & seq(i+1)=='U'
            count(26)=count(26)+1;
        elseif b(i:i+2)=='.(.' & seq(i+1)=='C'
            count(27)=count(27)+1;
        elseif b(i:i+2)=='.(.' & seq(i+1)=='G'
            count(28)=count(28)+1;
        elseif b(i:i+2)=='.((' & seq(i+1)=='A'
            count(29)=count(29)+1;
        elseif b(i:i+2)=='.((' & seq(i+1)=='U'
            count(30)=count(30)+1;
        elseif b(i:i+2)=='.((' & seq(i+1)=='C'
            count(31)=count(31)+1;
        elseif b(i:i+2)=='.((' & seq(i+1)=='G'
            count(32)=count(32)+1;
        end
        feature(k,11:42)=count;
    end
    
end
for i=1:n
    base=basecount(cell2mat(SS(i,2)));
    base=cell2mat(struct2cell(base));
    feature(i,43:46)=base(1:4,:);
    dimer=dimercount(cell2mat(SS(i,2)));
    dimer=cell2mat(struct2cell(dimer));
    feature(i,47:62)=dimer(1:16,:);
    codon=codoncount(cell2mat(SS(i,2)));
    codon=cell2mat(struct2cell(codon));
    feature(i,63:126)=codon(1:64,:);
end