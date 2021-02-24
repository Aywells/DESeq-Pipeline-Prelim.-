raw_table = readtable('Raw_cl.csv');
logic = isnan(table2array(raw_table(:,3:15)));
[num,~] = size(raw_table);

for i = 1:10842
    for j = 1:13
        
        hold_view = logic(i,j);
        if hold_view == true
            raw_table(i,j+2) = (array2table(0));
        end
        
    end
end

disease_pheno = raw_table(:,3:15);
gene_count = raw_table(:,19:4653);
patient_archive = cell(1,4653);

new_dis_pheno = [];

for i = 1:num
    for j = 1:13
        
        if table2array(disease_pheno(i,j)) == 1
            new_dis_pheno(i,:) = j;
        elseif sum(table2array(disease_pheno(i,:))) == 0
            new_dis_pheno(i,:) = 13;
        end
        
    end
end

gene_count_new = normalize(gene_count);
    [coeff,score,latent] = pca((table2array(gene_count_new))); % Load PCA -(after transpose): row=artwork,column=PCs
    tPCAscore = score(:,1:3);

 %   scatter3(tPCAscore(:,1),tPCAscore(:,2),tPCAscore(:,3))

combine = [new_dis_pheno,tPCAscore(:,1:2)];
    
cmap = colormap(hsv(13));
hold on
for k1 = 1:num
    scatter(tPCAscore(k1,1),tPCAscore(k1,2),[4],cmap(new_dis_pheno(k1),:))
end
[tPCAscore(:,1),tPCAscore(:,2),new_dis_pheno];
hold on
k1 = 1
scatter(tPCAscore(k1,1),tPCAscore(k1,2),5,cmap(new_dis_pheno(3),:))
k1 = 2
scatter(tPCAscore(k1,1),tPCAscore(k1,2),5,cmap(new_dis_pheno(3),:))
k1 = 3
scatter(tPCAscore(k1,1),tPCAscore(k1,2),5,cmap(new_dis_pheno(13),:))




[idx,C] = kmeans(table2array(tPCAscore(:,1:2)),13);
hold on
plot(tPCAscore(:,1),tPCAscore(:,2),'b.','MarkerSize',5)
plot(C(:,1),C(:,2),'r.','MarkerSize',20)


lda = FDA((table2array(gene_count))',new_dis_pheno,2);
plot(lda(1,:)',lda(2,:)','.')
cmap = colormap(hsv(13));
hold on
for k1 = 1:num
    scatter(lda(1,:)',lda(2,:)',[4],cmap(new_dis_pheno(k1),:))
end


svm = fitcecoc(gene_count, new_dis_pheno);