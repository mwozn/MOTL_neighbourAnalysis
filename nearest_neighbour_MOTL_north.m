% Measure distance to nearest neighbour from MOTL input's centre
% Use subtom_transform_motl.sh to re-centre
% Michael Wozny 20210412

% path to allmotl file to start from:
input_motl_fn = 'allmotl_clean_north_16.em';

% read the input motl to start from
input_motl = getfield(tom_emread(input_motl_fn), 'Value');

% tomogram number, MOTL row 6
tomogram_number = unique(input_motl(6,:));

for tom_num = 1:length(tomogram_number)
    % tomogram_subtomos holds information of subtomos according to row 6
    tomogram_subtomos = [];
    for k = 1:length(input_motl)
        if ismember(input_motl(6,k),tomogram_number(tom_num))
            subtomo = input_motl(:,k);
            if isempty(subtomo) == 1
                continue
            elseif isempty(tomogram_subtomos) == 1
                tomogram_subtomos = subtomo;
                continue
            else
                tomogram_subtomos = horzcat(tomogram_subtomos,subtomo);
            end
        end
    end
    
    % subtomo_centers holds xyz coordinates after translation of
    % tomogram_subtomos
    [r,q] = size(tomogram_subtomos);
    subtomo_centers = [];
    for k = 1:q
        % row 1 = x after translation
        subtomo_centers(1,k) = tomogram_subtomos(8,k) + tomogram_subtomos(11,k);
        % row 2 = y after translation
        subtomo_centers(2,k) = tomogram_subtomos(9,k) + tomogram_subtomos(12,k);
        % row 3 = z after translation
        subtomo_centers(3,k) = tomogram_subtomos(10,k) + tomogram_subtomos(13,k);
    end
    
    % loop over model_table and get the Euclidean norm difference between each
    % point, store these in norm_dist and then find the minimum distance for each
    % point
    norm_dist = [];
    [r,q] = size(subtomo_centers);
    for k = 1:q
        if size(subtomo_centers,2) == 2
            norm_dist = norm(subtomo_centers(:,1) - subtomo_centers(:,2));
            break
        end
        for kk = 1:q
            norm_dist(k,kk) = norm(subtomo_centers(:,k) - subtomo_centers(:,kk));
        end
    end
    
    % collect the nearest neighbour of each dipole for each model, in px
    nearest_neighbour = [];
    for k = 1:length(norm_dist)
        subtomo_norm_dist = norm_dist(:,k);
        if size(norm_dist,1) == 1
            break
        end
        nearest_neighbour(k) = min(subtomo_norm_dist(subtomo_norm_dist>0));
    end
    
    % round to 4 digits and transpose arrary
    nearest_neighbour = transpose(nearest_neighbour);
    nearest_neighbour = round(nearest_neighbour,4);

    % report in nm
    pixel_size = 0.2684;
    bin = 2;
    NEAREST_NEIGHBOUR{tom_num} = (nearest_neighbour * (pixel_size * bin));
end

% concatenate NEAREST_NEIGHBOUR_ALL to have all nearest neighbour
% measurements together
all = cat(1, NEAREST_NEIGHBOUR{:});
disp(strcat(['Median = ',num2str(median(all))]));

% report bin with max counts
binWidth = 1;
counts = histcounts(all,'BinWidth',binWidth,'BinLimits',[0,100]);
[row,col]=find(ismember(counts,max(counts)));
disp(strcat(['Hist peak = ',int2str((col * binWidth)-binWidth),'-',int2str(col * binWidth)]));

csvwrite('nearest_neighbour_MOTL_north.txt',all)

% plot nearest neighbour measurements as a histogram
figure('Renderer', 'painters', 'Position', [10 10 1200 600])
histogram(all,'BinWidth',binWidth,'BinLimits',[0,100])
xlabel('Nearest Dipole Center (nm)')
set(gca,'XTick',0:5:100)
ylim([0 45])