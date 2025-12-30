y_ctd = all_data.year(all_data.type==1);
y_osd = all_data.year(all_data.type==2);
y_flt = all_data.year(all_data.type==3);

counts_ctd = histc(y_ctd,1965:2025);
counts_osd = histc(y_osd,1965:2025);
counts_flt = histc(y_flt,1965:2025);

figure;
bar(1965:2025,[counts_ctd counts_osd counts_flt],'stacked');
grid
