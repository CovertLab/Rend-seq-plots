%% import .wig and compute z score
% Import data and compute peak and step z score.

% parameters for peak z score and step z score calculation
half_width_z_score = 100;    % half width of window for averaging for peak z score
half_width_step = 100;      % half width of window for averaging for step z score
average_threshold = 0.25;   % average read density (read/nt) threshold for consideration
gap_z = 2;                  % gap left out (both sides) of central position for peak z score 
gap_step = 3;               % gap left out (both sides) of central position for step z score 

cache = 0;
cache_file = 'cache.mat';
cache_vars = {'data', 'z_peak', 'z_step'};

% full path of the directory containing the wig files. Needs to be modified
% to directory of interest
data_dir = 'data/wigs_U000962';
addpath('sub_routines')

%% 5' exo WT read file
files{1} = 'GSM2971252_Escherichia_coli_WT_Rend_seq_5_exo_MOPS_comp_25s_frag_pooled_*_no_shadow.wig';

% genome size in nucleotides:
genome_size = 4639675;

% get_data_compute_statistic_simple: imports the wig file in a Matlab array
% and computes the peak and z score. Takes about ~30 s per dataset.
if cache
	load(cache_file, cache_vars{:});
else
	[data,z_peak,z_step] =  get_data_compute_statistic(...
	    files,data_dir,...
	    half_width_z_score,half_width_step,...
	    average_threshold,gap_step,gap_z,genome_size);
	save(cache_file, cache_vars{:});
end
%%
%% Remove this part later:
annotation_file = 'U00096.2.faa';
annotation_dir = 'genome';
[start_f, stop_f, genes_f, start_r, stop_r, genes_r] = read_gene_annotation_20180206(annotation_file,annotation_dir);

window_size = 10;
frac_reads = 0.5;
gene_pad = 100;

% get start and stop of regions split by segments of no reads for fwd strand
fwd_reads = data(1,:,1) + data(1,:,3);
points = 1:length(fwd_reads)-window_size;
reads = zeros(length(points), 1);
for i = points
	reads(i) = sum(fwd_reads(i:i+window_size) > 0) > frac_reads * window_size;
end

real_starts = [];
real_stops = [];
gene = 1;
pos = 1;

shifts = points(reads == 0);
n_genes = length(start_f);
n_shifts = length(shifts);
while gene <= n_genes
	if gene == 1 & start_f(gene) < shifts(pos)
		real_starts = [real_starts 1]
	else
		real_starts = [real_starts shifts(pos)];
	end

	while gene < n_genes & pos < n_shifts
		if stop_f(gene) > shifts(pos) + window_size
			pos = pos + 1;
		elseif shifts(pos) > start_f(gene + 1) - gene_pad
			gene = gene + 1;
		else
			temp_pos = pos;
			next_start = start_f(gene + 1);
			while temp_pos < n_shifts & shifts(temp_pos + 1) < next_start
				temp_pos = temp_pos + 1;
			end
			real_stops = [real_stops shifts(temp_pos)];
			break
		end
	end
	gene = gene + 1;
end

% last shift is before end of last gene so set end at end of genome
if length(real_starts) ~= length(real_stops)
	real_stops = [real_stops genome_size];
end
%%

%Remove this section as well: Just here for testing right now:

strand_region = 1;
a=b
for i = 1:length(real_starts)
	i
	start_region = real_starts(i);
	stop_region = real_stops(i);

	% identification, quantification and plot of region of interest.
	[x_3, x_5, isoform_levels] = quantify_plot_isoforms_20180130(start_region,stop_region,...
	    strand_region,z_thresh,width_z,z_peak_oi,data_oi,known_5,known_3,...
	    spurious_5,spurious_3,genome_size,annotation_file,annotation_dir,strcat('wt/fwd/', num2str(i)));
	% outputs of above:
	% x_3: position of 3' ends in region of interest.
	% x_5: position of 5' ends in region of interest.
	% isoform_levels: isoform level array, entry at (i,j) equals abundance of
	%                 mRNA isoform starting at x_5(i) and ending at x_3(j).

	if length(x_3) > 0 & length(x_5) > 0
		operon_start = min(x_5);
		operon_stop = max(x_3);
		gene_index = (start_f > operon_start) & (stop_f < operon_stop);
		genes_in_isoforms = genes_f(gene_index);
		gene_start = start_f(gene_index);
		gene_stop = stop_f(gene_index);
		n_genes = length(genes_in_isoforms);

		added_tu = false;
		for j = 1:length(x_3)
			for k = 1:length(x_5)
				level = isoform_levels(k, j);
				if level > 0
					isoform_start = x_5(k);
					isoform_stop = x_3(j);
					genes = genes_in_isoforms(gene_start > isoform_start & gene_stop < isoform_stop);
					if length(genes) > 0
						gene_str = genes{1};
						for m = 2:length(genes)
							gene_str = strcat(gene_str, ', ', genes{m});
						end
						fprintf(output_file, '%d\t%d\t"%s"\t%d\t%d\t%d\t%0.4f\t"%s"\n', operon, n_genes, 'f', i, isoform_start, isoform_stop, level, gene_str);
						added_tu = true;
					end
				end
			end
		end

		if added_tu
			operon = operon + 1;
		end
	end
end


% get start and stop of regions split by segments of no reads for rev strand
new_start_r = sort(genome_size - stop_r);
new_stop_r = sort(genome_size - start_r);

rev_reads = data(1,:,2) + data(1,:,4);
points = 1:length(rev_reads)-window_size;
reads = zeros(length(points), 1);
for i = points
	reads(i) = sum(rev_reads(i:i+window_size) > 0) > frac_reads * window_size;
end

real_starts = [];
real_stops = [];
gene = 1;
pos = 1;

shifts = points(reads == 0);
n_genes = length(new_start_r);
n_shifts = length(shifts);
while gene <= n_genes
	if gene == 1 & new_start_r(gene) < shifts(pos)
		real_starts = [real_starts 1]
	else
		real_starts = [real_starts shifts(pos)];
	end

	while gene < n_genes & pos < n_shifts
		if new_stop_r(gene) > shifts(pos) + window_size
			pos = pos + 1;
		elseif shifts(pos) > new_start_r(gene + 1) - gene_pad
			gene = gene + 1;
		else
			temp_pos = pos;
			next_start = new_start_r(gene + 1);
			while temp_pos < n_shifts & shifts(temp_pos + 1) < next_start
				temp_pos = temp_pos + 1;
			end
			real_stops = [real_stops shifts(temp_pos)];
			break
		end
	end
	gene = gene + 1;


end

% last shift is before end of last gene so set end at end of genome
if length(real_starts) ~= length(real_stops)
	real_stops = [real_stops genome_size];
end
