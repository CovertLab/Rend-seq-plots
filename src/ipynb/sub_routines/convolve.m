function filtered_data = convolve(data,filter)

% discrete Fourier transform of the data and of the filter.
ft_data = fft(data);
ft_filter = fft(filter);

    filtered_data = fftshift(ifft(ft_filter'.*ft_data));
end