function [data_FinalRules,YR_table]=FasTApriori(Numcore,T,T_X,T_Y,MinSupp1,MinSupp2,MinConf)
%% Copyright (c) 2021, Linwei Li
%% Project Title: An improved Apriori algorithm for deformation response analysis of landslides
%% Developer: Linwei Li
%% Contact Info: linweili93@126.com, lwli@gzu.edu.cn
% T={T_X,T_Y}
%% Obtain the the front itemset IF from the transaction database T
    IF_Items=[];
    for i=1:numel(T_X)
        IF_Items=union(IF_Items,T_X{i});
    end
    IF_Items=IF_Items(:)';
%% Obtain the the rear itemset IR from the transaction database T
    IR_Items=[];
    for i=1:numel(T_Y)
        IR_Items=union(IR_Items,T_Y{i});
    end
    IR_Items=IR_Items(:)';      
%% the 1st scanning of the transaction database to filter the invalid items in IF
    Temp_IF=cell(1,numel(IF_Items));
    IF_count=num2cell(zeros(1,size(Temp_IF,2)));
    for r=1:numel(IF_Items)
        Temp_IF{1,r}=IF_Items{1,r};
        for i=1:numel(T_X)
            if all(ismember(Temp_IF{1,r},T_X{i}))
                IF_count{1,r}=IF_count{1,r}+1;
            end
        end
    end
    IF_Items=Temp_IF(cell2mat(IF_count)/numel(T_Y)>=MinSupp1);
%% the 1st scanning of the transaction database to filter the invalid items in IR  
    Temp_IR=cell(1,numel(IR_Items));
    IR_count=num2cell(zeros(1,size(Temp_IR,2)));
    for r=1:numel(IR_Items)
        Temp_IR{1,r}=IR_Items{1,r};
        for i=1:numel(T_Y)
            if all(ismember(Temp_IR{1,r},T_Y{i}))
                IR_count{1,r}=IR_count{1,r}+1;
            end
        end
    end
    IR_Items=Temp_IR(cell2mat(IR_count)/numel(T_Y)>=MinSupp1);
%% The eligible items were stored in the effective front itemset IF_S by different dimensions and levels.
    i=1;
    while i<=length(IF_Items)
        if i==1
        temp_char=abs(char(IF_Items{i}));
        elseif sum(Temp_IF_level)<length(IF_Items)
        temp_char=abs(char(IF_Items{sum(Temp_IF_level)+1}));
        elseif sum(Temp_IF_level)>=length(IF_Items)
        break
        end
        if temp_char(2)==48
            temp_name=strcat('F0',strcat(num2str(i),'-'));
        else
            temp_name=strcat('F',strcat(num2str(i),'-'));
        end
        if sum(cell2mat(strfind(IF_Items,temp_name)))~=0
            Temp_IF_level(1,i)=sum(cell2mat(strfind(IF_Items,temp_name)));
            Temp_IF_index=~cellfun('isempty',strfind(IF_Items,temp_name));
            IF_Dim_level{1,i}=IF_Items(Temp_IF_index);
            i=i+1;
        else
            break    
        end
    end
    IF_S=IF_Dim_level;
%% The eligible items were stored in the effective rear itemset IR_S by different dimensions and levels.
    i=1;
    while i<=length(IR_Items)
        if i==1
            temp_char=abs(char(IR_Items{i}));
        elseif sum(Temp_IR_level)<length(IR_Items)
            temp_char=abs(char(IR_Items{sum(Temp_IR_level)+1}));
        elseif sum(Temp_IR_level)>=length(IR_Items)
            break
        end
        if temp_char(2)==48
            temp_name=strcat('R0',strcat(num2str(i),'-'));
        else
            temp_name=strcat('R',strcat(num2str(i),'-'));
        end
        if sum(cell2mat(strfind(IR_Items,temp_name)))~=0
            Temp_IR_level(1,i)=sum(cell2mat(strfind(IR_Items,temp_name)));
            Temp_IR_index=~cellfun('isempty',strfind(IR_Items,temp_name));
            IR_Dim_level{1,i}=IR_Items(Temp_IR_index);
            i=i+1;
        else
            break
        end
    end
    IR_S=IR_Dim_level;
%%  Generate the antecedent set XF through the iterative computation of Cartesian products between different dimensions of items stored in IF_S
    All_Comb=ff2n(size(IF_S,2));
    All_Comb(1,:)=[];
    temp_All_Items=cell(size(All_Comb,1),1);
    
    parfor i=1:size(All_Comb,1)
    temp_IF_S=IF_S(logical(All_Comb(i,:)));
    if size(temp_IF_S,2)~=1
    k=1;
    while k<=size(temp_IF_S,2)-1
        if k==1
            [x, y]=meshgrid(temp_IF_S{1,k},temp_IF_S{1,k+1});
            temp_X=reshape(x,[1,size(x,1)*size(x,2)]);
            temp_Y=reshape(y,[1,size(y,1)*size(y,2)]);
            temp_Items=arrayfun(@(k) {temp_X{k},temp_Y{k}},1:length(temp_Y),'Un',0)';
        else
            [x, y]=meshgrid(temp_Items,temp_IF_S{1,k+1});
            temp_X=reshape(x,[1,size(x,1)*size(x,2)]);
            temp_Y=reshape(y,[1,size(y,1)*size(y,2)]);
            temp_Items=arrayfun(@(a){temp_X{a}{1:k},temp_Y{a}},1:length(temp_Y),'Un',0)';
        end
        k=k+1;
    end
    else
       temp_Items=temp_IF_S{:}';
    end
       temp_All_Items{i,1}=temp_Items;
    end
    
    k=1;
    for i=1:size(temp_All_Items,1)
        for j=1:size(temp_All_Items{i,1},1)
            XF{k,1} = {temp_All_Items{i}{j,1}};
            k=k+1;
        end
    end
%%  Generate the antecedent set YR through the iterative computation of Cartesian products between different dimensions of items stored in IR_S
    k=1;
    while k<=size(IR_S,2)
        if k==1
            temp_R_Items=(IR_S{:})';
        else
            [x, y]=meshgrid(temp_R_Items,IR_S{k+1});
            temp_X=reshape(x,[1,size(x,1)*size(x,2)]);
            temp_Y=reshape(y,[1,size(y,1)*size(y,2)]);
            temp_R_Items=arrayfun(@(a) {temp_X{a}{1:k},temp_Y{a}},1:length(temp_Y),'Un',0)';
        end
        k=k+1;
    end
    YR=temp_R_Items;
%%  the 2nd scanning of the transaction database to otain the effective antecedent set XF_S  
    Multiplier_matrix=10.^(size(XF{1}{1},2)-1:-1:0);
    order=1:length(T_X);
    order1=buffer(order,Numcore);
    spmd
    Temp_XF_count_matrix=sparse(size(XF,1),sum(order1(labindex,:)~=0));
    for i=1:sum(order1(labindex,:)~=0)
        T_ASCII=abs(char(T_X{order1(labindex,i)}(:)));
        T_encode=sum(Multiplier_matrix.*T_ASCII,2);
        for k=1:size(XF,1)
                XF_ASCII=abs(char(XF{k}{:}));
                XF_encode=sum(Multiplier_matrix.*XF_ASCII,2);
                T_encode_matrix= (T_encode*ones(1,size(XF_ASCII,1)))';
                Temp_XF_count_matrix(k,i)=double(sum(sum(T_encode_matrix-XF_encode==0))==size(XF_ASCII,1));             
        end
%         1+i/numel(T_X)
    end
    end
    XF_count_matrix=sparse(size(XF,1),size(T_X,1));
    for i=1:size(order1,1)
        XF_count_matrix(:,order1(i,order1(i,:)~=0))=Temp_XF_count_matrix{1,i};
    end
    Temp_XF_count_matrix=XF_count_matrix;
    XF_total_count=sum((Temp_XF_count_matrix),2);
    XF_S=XF(XF_total_count/numel(T_X)>=MinSupp1);
    XF_S_count=XF_total_count(XF_total_count/numel(T_X)>=MinSupp1); 
%%  the 2nd scanning of the transaction database to otain the effective antecedent set YR_S  
% Since the number of the rear items is generally significantly less than 
% that of the front items, different calculation method was adopted here
    Temp_YR_count=zeros(length(YR),1); 
    for i=1:numel(T)
        for r=1:numel(YR)
                Temp_YR_count(r)=Temp_YR_count(r)+1*double(size(unique(cat(2,T{i},YR{r})),2)==size(char(T{i}),1));
        end
    end
    YR_total_count=Temp_YR_count;
    YR_S=YR(YR_total_count/numel(T_Y)>=MinSupp1);
    YR_S_count=YR_total_count(YR_total_count/numel(T_Y)>=MinSupp1); 
    YR_rank=1:length(YR_S);
    YR_table=cell2table(YR_S);
    YR_table.Rcount=YR_S_count;
%% Generate the potential transaction set V
    [x, y]=meshgrid(XF_S,YR_S);
    X=reshape(x,[1,length(x)*size(x,1)]);
    Y=reshape(y,[1,length(y)*size(y,1)]);
    V_Items=arrayfun(@(k) [x{k},y{k}],1:length(x)*size(x,1),'Un',0)';
    [~, temp_rank]=meshgrid(XF_S,YR_rank);
    Y_rank=reshape(temp_rank,[1,length(temp_rank)*size(temp_rank,1)]);
%%   the 3rd scanning of the transaction database to obtain the strong association rules (SARule)
    order=1:length(T);
    order2=buffer(order,Numcore);
    spmd
    Temp_V_count_matrix=sparse(numel(V_Items),sum(order2(labindex,:)~=0));
    for i=1:sum(order2(labindex,:)~=0)
        T_ASCII=abs(char(T{order2(labindex,i)}{:}));
        T_encode=sum(Multiplier_matrix.*T_ASCII,2);
        for k=1:size(V_Items,1)
            Temp_V={V_Items{k}{:}};
            Temp_Vitems={};
            o=1;
            for j=1:size(Temp_V,2)
                if iscell(Temp_V{1,j})
                    for p=1:size(Temp_V{1,j},2)
                        Temp_Vitems{1,o}=Temp_V{1,j}{p};
                        o=o+1;
                    end
                else
                    Temp_Vitems{1,o}=Temp_V{1,j};
                    o=o+1;
                end
            end
            V_ASCII=abs(char(Temp_Vitems)); 
            V_encode=sum(Multiplier_matrix.*V_ASCII,2);
            T_encode_matrix= (T_encode*ones(1,size(V_ASCII,1)))'; 
            Temp_V_count_matrix(k,i)=double(sum(sum(T_encode_matrix-V_encode==0))==size(V_ASCII,1)); 
         end
%         2+i/numel(T)
    end
    end
    v_count_matrix=sparse(numel(V_Items),size(T,1));
    for i=1:size(order2,1)
        v_count_matrix(:,order2(i,order2(i,:)~=0))=Temp_V_count_matrix{1,i};
    end
    Temp_V_count_matrix=v_count_matrix;

    V_total_count=sum((Temp_V_count_matrix),2);
    V_tureItems=V_Items(V_total_count/numel(T)>=MinSupp2);
    V_count=V_total_count(V_total_count/numel(T)>=MinSupp2);
    V_Supp=V_count/numel(T);
  
    SARule_Y=Y(V_total_count/numel(T)>=MinSupp2);
    SARule_Y=SARule_Y';
    SARule_X=X(V_total_count/numel(T)>=MinSupp2);
    SARule_X=SARule_X';
    SARule_Y_rank=Y_rank(V_total_count/numel(T)>=MinSupp2);
    SARule_Y_rank=SARule_Y_rank';

    [temp_x, temp_y]=meshgrid(XF_S_count,YR_S_count);
    Rule_X_count=reshape(temp_x,[1,length(temp_x)*size(temp_x,1)]);
    Rule_X_count=Rule_X_count(V_total_count/numel(T)>=MinSupp2)';
    Rule_X_Supp=Rule_X_count/numel(T);
    Rule_Y_count=reshape(temp_y,[1,length(temp_y)*size(temp_y,1)]);
    Rule_Y_count=Rule_Y_count(V_total_count/numel(T)>=MinSupp2)';
    Rule_Y_Supp=Rule_Y_count/numel(T);

    V_temp_conf=V_Supp./Rule_X_Supp;
    Rule_X2Y_X=SARule_X(V_temp_conf>=MinConf);
    Rule_X2Y_Y=SARule_Y(V_temp_conf>=MinConf);
    Rule_X2Y_V=V_tureItems(V_temp_conf>=MinConf);
    Rule_X2Y_Yrank=SARule_Y_rank(V_temp_conf>=MinConf);
    Rule_X2Y_V_Supp=V_Supp(V_temp_conf>=MinConf);
    Rule_X2Y_XSupp=Rule_X_Supp(V_temp_conf>=MinConf);
    Rule_X2Y_YSupp=Rule_Y_Supp(V_temp_conf>=MinConf);
    Rule_X2Y_Conf=V_temp_conf(V_temp_conf>=MinConf);
%% Sort According to Confidence
    FinalRules=cell(length(Rule_X2Y_Conf),5);
    temp_rank=cellfun('length',Rule_X2Y_V)-size(IR_S,2);
    temp_diff=max(temp_rank)-temp_rank;
    temp_order=temp_diff==max(temp_rank)-1;
    
    Rule_X2Y_X(1:sum(temp_order))=arrayfun(@(a) cat(2,Rule_X2Y_X{a},cell(1,temp_diff(a))),1:sum(temp_order),'Un',0)';
    Rule_X2Y_X(sum(temp_order)+1:length(Rule_X2Y_X))=arrayfun(@(a) cat(2,Rule_X2Y_X{a}{1:temp_rank(a)},cell(1,temp_diff(a))),sum(temp_order)+1:length(Rule_X2Y_X),'Un',0)';
    
    FinalRules(:,1)=Rule_X2Y_X;
    FinalRules(:,2)=Rule_X2Y_Y;
    FinalRules(:,3)=num2cell(Rule_X2Y_Yrank);
    FinalRules(:,4)=num2cell(Rule_X2Y_V_Supp);
    FinalRules(:,5)=num2cell(Rule_X2Y_Conf);

    data_FinalRules=cell(length(IR_Items),1);
    for i=1:length(IR_Items)
        temp_data=FinalRules(cell2mat(FinalRules(:,3))==i,:);
        [~, temp_SortOrder]=sort(cell2mat(temp_data(:,5)),'descend');
        temp_data2=temp_data(temp_SortOrder,:);
        data_FinalRules{i}=temp_data2;
    end
end