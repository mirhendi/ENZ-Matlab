%% Salam
clc;clear all;close all;
warning off;
%% Sequances
[Header ,Seq] = fastaread('input_sequences.txt');
%% Enzyme Database
fid = fopen('rebase_e.txt','r');
data = textscan(fid,'%s%s%f%f%f%f%f%f%f','headerlines',42,'delimiter','\t');
[Enzymes, ~, ~, ~,~,~,~,~,~] = deal(data{:});
fclose(fid);

%% Restriction Enzymes
result(75,numel(Seq),numel(Enzymes))=0;
for i=[1:684,687:numel(Enzymes)]
    for j=1:numel(Seq)
        [~, ~, Lengths]  = restrict(Seq{j},Enzymes{i});
        Lengths(Lengths<50)=[]; %Remove Small Seq.
        Lengths=round(10*log(Lengths));
        result(Lengths,j,i)=1;
    end
end


for i=1:numel(Enzymes)
    all_var=0;
    for j=1:size(result,1)
        all_var=all_var+var(result(j,:,i));
    end
    enz_var(i)=all_var;
end
[enz_var,enz_index]=sort(enz_var, 'descend' );

%% Marker
Marker=zeros(size(result,1),1);
Marker(round(10*log([100 200 300 400 500 600 700 800 900 1000 1500])))=1;


%% Gel
n=5; % number of best enzymes to shown
for i=1:n
    pic=draw_pic([Marker,result(:,:,enz_index(i))]);
    figure;imshow(pic,[]);title(Enzymes{enz_index(i)});
    saveas(gcf,['no_' num2str(i) '_' Enzymes{enz_index(i)}], 'jpg')
    disp(Enzymes{enz_index(i)});
    for j=1:numel(Seq)
    [~, ~, Lengths]=restrict(Seq{j},Enzymes{enz_index(i)});
    end
end

%% Make Table
enz_per_seq{numel(Seq),n}=0;
for i=1:n
    res_len=[];
    for j=1:numel(Seq)
        [~, ~, Lengths]=restrict(Seq{j},Enzymes{enz_index(i)});
        Lengths=sort(Lengths');
        res_len=[res_len;{int2str(Lengths)}];
    end
    enz_per_seq(:,i)=res_len;
end
Table = cell2table(enz_per_seq,'VariableNames', genvarname(Enzymes(enz_index(1:n))));
Table.Properties.RowNames = Header;
Table
% uf = uifigure;
% uitable(uf,'ColumnName',Header,'RowName',Enzymes(enz_index(1:n)),'Data',enz_per_seq,'ColumnWidth','auto')

