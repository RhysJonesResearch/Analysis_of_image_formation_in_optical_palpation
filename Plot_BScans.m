%%% Temp commands
%% Import OCT
Im1 = loadSSOCT(FilePath_Uncompressed,ImSize);
%%
Im2 = loadSSOCT(FilePath_Compressed,ImSize);

%% Plot
figure; 
subplot(2,8,1); imagesc(squeeze(Im1(:,1,:))'); colormap(gray); caxis([0 100]); title('y=1');
subplot(2,8,2); imagesc(squeeze(Im1(:,200,:))'); colormap(gray); caxis([0 100]); title('y=200');
subplot(2,8,3); imagesc(squeeze(Im1(:,400,:))'); colormap(gray); caxis([0 100]); title('y=400');
subplot(2,8,4); imagesc(squeeze(Im1(:,600,:))'); colormap(gray); caxis([0 100]); title('y=600');
subplot(2,8,5); imagesc(squeeze(Im1(:,800,:))'); colormap(gray); caxis([0 100]); title('y=800');
subplot(2,8,6); imagesc(squeeze(Im1(:,1000,:))'); colormap(gray); caxis([0 100]); title('y=1000');
subplot(2,8,7); imagesc(squeeze(Im1(:,1200,:))'); colormap(gray); caxis([0 100]); title('y=1200');
%subplot(2,8,8); imagesc(squeeze(Im1(:,1400,:))'); colormap(gray); caxis([0 100]); title('y=1400');
subplot(2,8,9); imagesc(squeeze(Im2(:,1,:))'); colormap(gray); caxis([0 100]); title('y=1');
subplot(2,8,10); imagesc(squeeze(Im2(:,200,:))'); colormap(gray); caxis([0 100]); title('y=200');
subplot(2,8,11); imagesc(squeeze(Im2(:,400,:))'); colormap(gray); caxis([0 100]); title('y=400');
subplot(2,8,12); imagesc(squeeze(Im2(:,600,:))'); colormap(gray); caxis([0 100]); title('y=600');
subplot(2,8,13); imagesc(squeeze(Im2(:,800,:))'); colormap(gray); caxis([0 100]); title('y=800');
subplot(2,8,14); imagesc(squeeze(Im2(:,1000,:))'); colormap(gray); caxis([0 100]); title('y=1000');
subplot(2,8,15); imagesc(squeeze(Im2(:,1200,:))'); colormap(gray); caxis([0 100]); title('y=1200');
%subplot(2,8,16); imagesc(squeeze(Im2(:,1400,:))'); colormap(gray); caxis([0 100]); title('y=1400');
set(gcf, 'Position',  [0, 0, 2000, 1000]);

%% Global remove Line
if config.RemoveLine == 1
    img = squeeze(Im1(:,1,:));
    m = mean2(img);
    nAv = 500;
    for n = 1:nAv
        %img = img + squeeze(Im2(:,1+n*50,:));
        img = img + squeeze(Im1(:,1+n*floor(ImSize(2)/nAv),:));
    end
    img = img / (nAv+1);

    figure
    imagesc(img')
    caxis('auto');
    colorbar;
    colormap gray;

    img(img<40) = 0;
    img(:,1:350) = 0;
    img(:,390:end) = 0;
    img(img>0) = m;
    img = medfilt2(img,[3 3]);
    figure
    imagesc(img')
    caxis('auto');
    colorbar;
    colormap gray;
    lineremove = img';

    key = lineremove > 0;
    key1 = key;
    % figure
    % imagesc(key1)
    % caxis('auto');
    % colorbar;
    % colormap gray;

    for n = 1:size(lineremove,2)
        ind = find(key(:,n));
        move = max(ind)-min(ind)+1;
        key(ind,n) = 0;
        if rem(move,2) == 0
            key(ind(1:move/2)-move/2,n) = 1;
            key(ind((move/2+1):end)+move/2,n) = 1;
        else
            key(ind(1:ceil(move/2))-ceil(move/2),n) = 1;
            key(ind((ceil(move/2)+1):end)+floor(move/2),n) = 1;
        end
        %key(ind-move,n) = 1;
    end

    key2 = key;

    figure
    imagesc(key2-key1)
    caxis('auto');
    colorbar;
    colormap gray;
    
    img = squeeze(Im1(:,1,:))';
    img(key1) = img(key2);
    
    figure
    imagesc(img)
    caxis('auto');
    colorbar;
    colormap gray;
end

%% Local Remove Line 
% n = 1;
% if config.RemoveLine == 1
%     img = medfilt2(squeeze(Im1(:,n,:))',[5 5]);
%     nAv = 20;
%     for n = 1:nAv
%         %img = img + squeeze(Im2(:,1+n*50,:));
%         img = img + medfilt2(squeeze(Im1(:,1+n,:))',[5 5]);
%     end
%     img = img / (nAv+1);
% 
%     figure
%     imagesc(img)
%     caxis('auto');
%     colorbar;
%     colormap gray;
% 
%     img(img<57) = 0;
%     img(1:450,:) = 0;
%     img(565:end,:) = 0;
%     img(img>0) = 1;
%     %img = medfilt2(img,[1 3]);
%     img = bwareaopen(img,2000,4);
% %     img = medfilt2(img,[1 3]);
% %     img = medfilt2(img,[3 1]);
%     figure
%     imagesc(img)
%     caxis('auto');
%     colorbar;
%     colormap gray;
%     lineremove = img;
% 
%     key = lineremove > 0;
%     key1 = key;
%     % figure
%     % imagesc(key1)
%     % caxis('auto');
%     % colorbar;
%     % colormap gray;
% 
%     for n = 1:size(lineremove,2)
%         ind = find(key(:,n));
%         move = max(ind)-min(ind)+1;
%         key(ind,n) = 0;
%         if rem(move,2) == 0
%             key(ind(1:move/2)-move/2,n) = 1;
%             key(ind((move/2+1):end)+move/2,n) = 1;
%         else
%             key(ind(1:ceil(move/2))-ceil(move/2),n) = 1;
%             key(ind((ceil(move/2)+1):end)+floor(move/2),n) = 1;
%         end
%         %key(ind-move,n) = 1;
%     end
% 
%     key2 = key;
% 
%     figure
%     imagesc(key2-key1)
%     caxis('auto');
%     colorbar;
%     colormap gray;
%     
%     img = squeeze(Im1(:,1,:))';
%     
%     figure
%     imagesc(img)
%     caxis('auto');
%     colorbar;
%     colormap gray;
%     
%     img(key1) = img(key2);
%     
%     figure
%     imagesc(img)
%     caxis('auto');
%     colorbar;
%     colormap gray;
% end
%% Plot
figure; 
subplot(2,8,1); imagesc(squeeze(Im1(:,1,:))'); colormap(gray); caxis([0 100]); title('y=1');
subplot(2,8,2); imagesc(squeeze(Im1(:,400,:))'); colormap(gray); caxis([0 100]); title('y=200');
subplot(2,8,3); imagesc(squeeze(Im1(:,600,:))'); colormap(gray); caxis([0 100]); title('y=400');
subplot(2,8,4); imagesc(squeeze(Im1(:,600,:))'); colormap(gray); caxis([0 100]); title('y=600');
subplot(2,8,5); imagesc(squeeze(Im1(:,800,:))'); colormap(gray); caxis([0 100]); title('y=800');
subplot(2,8,6); imagesc(squeeze(Im1(:,1000,:))'); colormap(gray); caxis([0 100]); title('y=1000');
subplot(2,8,7); imagesc(squeeze(Im1(:,1200,:))'); colormap(gray); caxis([0 100]); title('y=1200');
%subplot(2,8,8); imagesc(squeeze(Im1(:,1400,:))'); colormap(gray); caxis([0 100]); title('y=1400');
subplot(2,8,9); imagesc(squeeze(Im2(:,420,:))'); colormap(gray); caxis([0 100]); title('y=1');
subplot(2,8,10); imagesc(squeeze(Im2(:,470,:))'); colormap(gray); caxis([0 100]); title('y=200');
subplot(2,8,11); imagesc(squeeze(Im2(:,520,:))'); colormap(gray); caxis([0 100]); title('y=400');
subplot(2,8,12); imagesc(squeeze(Im2(:,570,:))'); colormap(gray); caxis([0 100]); title('y=600');
subplot(2,8,13); imagesc(squeeze(Im2(:,620,:))'); colormap(gray); caxis([0 100]); title('y=800');
subplot(2,8,14); imagesc(squeeze(Im2(:,670,:))'); colormap(gray); caxis([0 100]); title('y=1000');
subplot(2,8,15); imagesc(squeeze(Im2(:,720,:))'); colormap(gray); caxis([0 100]); title('y=1200');
%subplot(2,8,16); imagesc(squeeze(Im2(:,1400,:))'); colormap(gray); caxis([0 100]); title('y=1400');
set(gcf, 'Position',  [0, 0, 2000, 1000]);