%% MTF PROJECT IKONO
clear all;
clc;
close all;

%% Step 1: access the path to the image
ParentFolderPath = 'C:\Users\santa\MTF Project\ProjectMTF\Aug01.2005.IK.mtfcOff';
ImagePath = dir(fullfile(ParentFolderPath, '*.tif'));
[image, raster] = readgeoraster(fullfile(ImagePath(1).folder, ImagePath(2).name));

% if your image is multi-band, uncomment this to select the first band
% image=image(:,:,1);

% show the full image
figure
imagesc(image), colormap("gray");

%% Step 2: MASKING THE IMAGE USING BINARY MASKING
% mask out invalid pixels (where image == 0)
binaryimage = ones(size(image), "logical");
binaryimage(image == 0) = 0;

%% Step 3: ROI MASK
% ROI geo-coordinates
X = [678204, 678243, 678243, 678204];
Y = [4906604, 4906604, 4906589, 4906589];

% convert geo-coords to pixel coords (map2pix is old but works)
[row, column] = map2pix(raster, X, Y);
row = round(row);
column = round(column);
binaryroiimage = poly2mask(column, row, size(image, 1), size(image, 2));

%% Step 4: first roi masking
% combine masks and get the final ROI
roi_image = logical(binaryimage .* binaryroiimage);
[row_find, column_find] = find(roi_image == 1);

% just using the first 13 rows found
row_find_subset = row_find(1:13, 1);
Extracted_ROI = roi_image(min(row_find_subset):max(row_find_subset), min(column_find):max(column_find));

% crop image to the ROI
figure
edge_spread_figure = image(min(row_find_subset):max(row_find_subset), min(column_find):max(column_find));
imagesc(edge_spread_figure), colormap("gray");

% plot each scan line
figure
for line = 1:size(Extracted_ROI, 1)
    plot(edge_spread_figure(line, :), 'linestyle', 'none', 'MarkerSize', 15, 'Marker', '.');
    hold on
end
title('Edge Spread Profile For 3 Scans');
xlabel('Pixel')
ylabel('Digital Number Value (Rows)')
grid on; grid minor;
set(gca, 'FontSize', 15, 'Fontweight', 'bold'); set(gcf, 'color', 'white');
hold off;

%% Step 5: Initialize Variables
% variables for spline fitting
num_rows = size(Extracted_ROI, 1);
num_cols = size(Extracted_ROI, 2);
roi = double(edge_spread_figure);
edge_midpoints = NaN(num_rows, 1);
x = 1:num_cols;
fine_x = linspace(1, num_cols, 1000); % finer grid for spline

%% Step 6: Spline Fitting and Inflection Point Detection
% find sub-pixel edge location for each row using a spline fit
figure;
hold on;
colors = lines(num_rows);
for row_idx = 1:num_rows
    y = roi(row_idx, :);
    
    % fit spline and get derivatives
    pp = spline(x, y);
    y_fine = ppval(pp, fine_x);
    pp_der1 = fnder(pp, 1);
    pp_der2 = fnder(pp, 2);
    y_der2 = ppval(pp_der2, fine_x);
    
    % find inflection points (where 2nd derivative is zero)
    zero_crossings = find(diff(sign(y_der2)));
    potential_inflections = [];
    for i = 1:length(zero_crossings)
        x_left = fine_x(zero_crossings(i));
        x_right = fine_x(zero_crossings(i) + 1);
        try
            root = fzero(@(z) ppval(pp_der2, z), [x_left, x_right]);
            potential_inflections(end + 1) = root;
        catch
            continue;
        end
    end
    
    % use the inflection point with the steepest slope
    if ~isempty(potential_inflections)
        if length(potential_inflections) > 1
            slopes = ppval(pp_der1, potential_inflections);
            [~, idx] = max(abs(slopes));
            midpoint = potential_inflections(idx);
        else
            midpoint = potential_inflections;
        end
        edge_midpoints(row_idx) = midpoint;
    end
    
    % plot the fit
    plot(fine_x, y_fine, 'Color', colors(row_idx, :), 'LineWidth', 1.5, 'DisplayName', sprintf('Row %d', row_idx));
end

xlabel('Pixel');
ylabel('Digital Number Value');
title('Spline Fit for Each Row');
legend('Location', 'best');
grid on;
hold off;

%% Step 7: Display Transition Points in Command Window
% display results
disp('Estimated Transition Points (Midpoints) for Each Row:');
for row_idx = 1:num_rows
    if ~isnan(edge_midpoints(row_idx))
        fprintf('Row %2d: %.4f pixels\n', row_idx, edge_midpoints(row_idx));
    else
        fprintf('Row %2d: Transition point not found.\n', row_idx);
    end
end
if ~isnan(edge_midpoints(row_idx))
    fprintf('Row %2d: %.4f pixels\n', row_idx, edge_midpoints(row_idx));
else
    fprintf('Row %2d: Transition point not found.\n', row_idx);
end

% define edge angle
e_y = edge_midpoints;
theta_deg = 6;
theta_rad = deg2rad(theta_deg);
cos_theta = cos(theta_rad);

%% Step 8: Calculate z(x, y) for Each Pixel
% project all pixels to a line perpendicular to the edge -> this creates the ESF
z = zeros(num_rows, num_cols);
x_indices = repmat(1:num_cols, num_rows, 1);
e_matrix = repmat(e_y, 1, num_cols);
z = (x_indices - e_matrix) * cos_theta;

%% Step 9: Plot Intensity vs. z(x, y)
% flatten data into vectors
z_vector = z(:);
intensity_vector = roi(:);
valid_indices = ~isnan(z_vector) & ~isnan(intensity_vector);
z_vector = z_vector(valid_indices);
intensity_vector = intensity_vector(valid_indices);

% plot the raw, oversampled ESF
figure;
scatter(z_vector, intensity_vector, 20, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('z(x, y) (pixels)');
ylabel('Digital Number');
title('Oversampled but Irregularly spaced ESF');
grid on;

% sort for the next step
[z_sorted, sort_idx] = sort(z_vector);
intensity_sorted = intensity_vector(sort_idx);

%% Step 10: Apply Modified Savitzky-Golay Filter for Filtering and Resampling
% filter the ESF with Savitzky-Golay
poly_order = 4;
window_size_z = 4;
oversampling_factor = 1;
z_min = min(z_sorted);
z_max = max(z_sorted);
z_uniform = z_min:(1 / oversampling_factor):z_max;
esf_filtered = zeros(size(z_uniform));

for i = 1:length(z_uniform)
    current_z = z_uniform(i);
    window_lower = current_z - (window_size_z / 2);
    window_upper = current_z + (window_size_z / 2);
    idx = find(z_sorted >= window_lower & z_sorted <= window_upper);
    if length(idx) >= (poly_order + 1)
        p = fit(z_sorted(idx), intensity_sorted(idx), 'poly4', 'Normalize', 'on', 'Robust', 'on');
        esf_filtered(i) = p(current_z);
    else
        esf_filtered(i) = NaN;
    end
end
esf_filtered = fillmissing(esf_filtered, 'linear');

%% Step 11: 
% hold off;
% figure;
% plot(z_uniform, esf_filtered, 'b-', 'LineWidth', 2);

%% Step 12: Display the Filtered ESF Values (Optional)
disp('Filtered and Uniformly Sampled ESF Values:');
disp([z_uniform', esf_filtered']);

%% Step 13: Identify Background and Foreground Regions
% define regions for SNR calculation
region_width = 10;
[~, edge_idx] = min(abs(z_uniform - 0));
background_start = 1;
background_end = edge_idx - region_width;
if background_end < background_start
    error('Background region is too close to the edge. Adjust region_width or verify ESF.');
end
foreground_start = edge_idx + region_width;
foreground_end = length(esf_filtered);
if foreground_start > foreground_end
    error('Foreground region is too close to the edge. Adjust region_width or verify ESF.');
end

% Extract background (dark) and foreground(bright) intensity values
background_intensity = esf_filtered(foreground_start:foreground_end);
foreground_intensity = esf_filtered(background_start:background_end);

%% Step 14: Calculate Edge Height
mean_background = mean(background_intensity);
mean_foreground = mean(foreground_intensity);
edge_height = mean_foreground - mean_background;

%% Step 15: Calculate Standard Deviations
std_background = std(background_intensity);
std_foreground = std(foreground_intensity);
avg_std = (std_background + std_foreground) / 2;

%% Step 16: Compute SNR
SNR = edge_height / avg_std;

% Plot filtered ESF over raw data
figure;
scatter(z_sorted, intensity_sorted, 20, 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
plot(z_uniform, esf_filtered, 'r-', 'LineWidth', 2);
xlabel('z(x, y) (pixels)');
ylabel('Digital Number');
title('Irregularly Spaced ESF and Filtered ESF');
legend('Original Data', 'Filtered ESF');
grid on;

% add SNR text box
annotation_text = sprintf('SNR: %.2f\nAverage Std Dev : %.2f\n', SNR, avg_std);
dim = [0.55 0.3 0.3 0.15];
annotation('textbox', dim, 'String', annotation_text, 'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 10);
hold off;

%% Step 17: Calculate the Line Spread Function (LSF)
% LSF is the derivative of the ESF
lsf_raw = diff(abs(esf_filtered));
lsf_raw(end + 1) = lsf_raw(end);
lsf_abs = abs(lsf_raw);
lsf_normalized = lsf_abs / max(lsf_abs);

% smooth LSF with a spline to get a better FWHM
z_fine = linspace(min(z_uniform), max(z_uniform), 10 * length(z_uniform));
lsf_spline = spline(z_uniform, lsf_normalized, z_fine);

% Plot LSF
figure;
plot(z_uniform, lsf_normalized, 'b-', 'LineWidth', 2);
xlabel('z(x, y) (pixels)');
ylabel('Normalized Absolute LSF');
title('Normalized Absolute Line Spread Function (LSF)');
grid on; grid minor;
hold on;

% trim LSF to a window around the peak for FWHM calc
hwidth = 100;
shift = 0;
[val, pos] = max(lsf_spline);
spoint = pos - hwidth + shift;
epoint = pos + hwidth + shift - 1;
spoint = max(spoint, 1);
epoint = min(epoint, length(z_fine));
t_real = lsf_spline(spoint:epoint);
t_real = t_real ./ max(t_real);
t_x = z_fine(spoint:epoint);

% find FWHM with linear interpolation
IHM = find(t_real > 0.5);
xl = (.5 - t_real(IHM(1))) * (-0.05) / (t_real(IHM(1) - 1) - t_real(IHM(1))) + t_x(IHM(1));
last = IHM(end);
xr = (.5 - t_real(last)) * (-0.05) / (t_real(last - 1) - t_real(last)) + t_x(last);

% plot FWHM lines
xline(xl, 'LineStyle', '--', 'Color', 'g', 'LineWidth', 1.5);
xline(xr, 'LineStyle', '--', 'Color', 'g', 'LineWidth', 1.5);
yline(0.5, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.5);

title('Line Spread Function (LSF) with FWHM');
ylabel('Normalized LSF Values (DN)');
xlabel('Pixel');
grid on; grid minor;
set(gca, 'FontSize', 15, 'FontWeight', 'bold');
set(gcf, 'color', 'white');
hold off;

%% Step 18: Calculate the Modulation Transfer Function (MTF) with Actual Spatial Frequencies
% MTF is the FFT of the LSF
mtf_fft = fft(lsf_normalized);
mtf_magnitude = abs(mtf_fft(1:floor(length(mtf_fft) / 2)));
mtf = mtf_magnitude / max(mtf_magnitude);

% define frequency axis
sampling_interval_pixels = 1 / oversampling_factor;
nyquist_freq_cpp = 1 / (2 * sampling_interval_pixels);
num_points = length(mtf);
frequency_axis_cpp = linspace(0, nyquist_freq_cpp, num_points);

% plot MTF
figure;
plot(frequency_axis_cpp, mtf, 'r-', 'LineWidth', 2);
xlabel('Spatial Frequency (cycles/pixel)');
ylabel('Normalized MTF');
title('Modulation Transfer Function (MTF) with Frequency in Cycles per Pixel');
xlim([0, nyquist_freq_cpp + 0.1]);
xticks(0:0.1:0.6)
ylim([0, 1.1]);
yticks(0:0.1:1.1);
grid on;
hold on;

% add reference lines
yline(0, '--k');
xline(nyquist_freq_cpp, '--k', sprintf('Nyquist Frequency = %.2f cycles/pixel', nyquist_freq_cpp), 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom');
mtf_at_nyquist = mtf(end);
plot(nyquist_freq_cpp, mtf_at_nyquist, 'bo', 'MarkerFaceColor', 'b');
text(nyquist_freq_cpp, mtf_at_nyquist, sprintf(' (%.3f, %.3f)', nyquist_freq_cpp, mtf_at_nyquist), 'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'right', 'Color', 'b');
hold off;