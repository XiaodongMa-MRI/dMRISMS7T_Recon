function sensmap=FitSensMap(sense_map,raw_data)

% created by Zhe Liu
% editted by Zhe Zhang

[nX,nY,nSL,nCH,nTE]=size(sense_map);

% get mask
SenMap.th = 0.80*1e-1; % Relative threshold for mask
% mask = 1.*(abs(ref)>SenMap.th.*max(abs(ref(:))));
for sli_idx = 1:nSL
    ref_img = sos(sqz(mean(raw_data(:,:,sli_idx,:,1),5)));
%     mask(:,:,sli_idx) = mmROI(ref_img,1,'roi',[0 SenMap.th.*max(abs(ref_img(:)))]);
    mask(:,:,sli_idx) = mmROI(ref_img,1,'roi',0);
end
% mask=get_mask(ref);
%     se = strel('disk',1);
%     mask=imdilate(imerode(mask,se),se);

%% Set up SenMap for sens. map assessment
SenMap.fit_method = 'Polynomial'; 
SenMap.full_size = [nX,nY];
% SenMap.Gauss_filter_kernel_width = 5;
SenMap.Gauss_filter_kernel_width = 5;
SenMap.Gauss_filter_kernel_sigma = 3;
SenMap.smooth_iter = 2;
% SenMap.fit_method = 'None';   
SenMap.fit_method = 'Polynomial' ;   
SenMap.fit_order = 6;
SenMap.fit_basis_num = 10;
SenMap.TPS_fit_method = 'indirect';
SenMap.ext_width = round((nX+nY)/100);          
SenMap.fit_margin = 0;
% SenMap.enlarge_fac = 2;


%% calc sens map
for sli_idx = 1:nSL
    nSL
    sens_raw = squeeze(sense_map(:,:,sli_idx,:));
    sensmap(:,:,:,sli_idx) = SensMapAssess(sens_raw, mask(:,:,sli_idx), SenMap);
end

end
