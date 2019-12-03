%%load(DiseaseSymScore.mat)
%load(miRNASecondStructure.mat)
%run(miRfeature.m)
z=size(Sym,1);
feature=miRfeature(SS);
for i=1:z
    a=max(Sym(i,:));
    Sym(i,:)=Sym(i,:)/a;    
end
[m,n]=find(MDA==1);
for i=1:length(m)
    P(i,1:126)=feature(m(i),:);
    P(i,127:448)=Sym(n(i),:);
end
pr=sum(P)/length(P);
[x,y]=find(MDA==0);
for i=1:length(x)
    UN(i,1:126)=feature(x(i),:);
    UN(i,127:448)=Sym(y(i),:);
end
for i=1:length(UN)
    d=0;
    for j=1:448
         d=d+(UN(i,j)-pr(j))^2;
         OD(i)=sqrt(d);
     end
end
averageOD=sum(OD)/length(UN);
z=find(OD>1.5*averageOD);
RN=UN(z,:);
RNsite=UNsite(z,:);
r=randperm(length(RN),length(P));
N=RN(r,:);
P(:,449)=1;
N(:,449)=0;
Sample=[P
    N];
INDICES = crossvalind('Kfold',length(P)*2,5);
for k=1:5
  test = (INDICES == k); 
  train = ~test;
  train_fea=Sample(train,1:448);
  train_labels=Sample(train,449);
  test_fea=Sample(test,1:448);
  test_target=Sample(test,449);
  mdl = TreeBagger(100, train_fea, train_labels);
  [pred_labels,score] = predict(mdl, test_fea);
  pred_labels=str2num(char(pred_labels'));
  C=confusionmat(test_target,pred_labels,'Order',[1 0]);
  TP=C(1,1);
  TN=C(2,2);
  FN=C(1,2);
  FP=C(2,1);
  acc(k)=(TP+TN)/sum(sum(C));
  sen(k)=TP/(TP+FN);
  spe(k)=TN/(TN+FP);
  pre(k)=TP/(TP+FP);
  mcc(k)=(TP*TN-FP*FN)/sqrt((TP+FN)*(TP+FP)*(TN+FN)*(TN+FP));
end
Acc=sum(acc)/5;
Spe=sum(spe)/5;
Sen=sum(sen)/5;
Pre=sum(pre)/5;
Mcc=sum(mcc)/5;
[X,Y,T,AUC]=perfcurve(test_target,score(:,2),1);
[XP,YP,TP,~] = perfcurve(test_target,score(:,2),1,'xcrit','reca','ycrit','prec');
YP(1)=YP(2);
AUCP=trapz(XP,YP);
