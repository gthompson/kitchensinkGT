function [w,scnl]=obspy_stream2waveform(path_to_python, path_to_converter, path_to_data)
%OBSPY_STREAM2WAVEFORM Use ObsPy to read a stream object and then convert
% into a waveform object. Also return a scnlobject.
% Example:
%  [w,scnl]=obspy_stream2waveform('/Users/glennthompson/anaconda/bin/python', ...
%       '/Users/glennthompson/obspy_stream2matfile.py', ...
%       'https://examples.obspy.org/BW.BGLD..EH.D.2010.037')
%  plot(w)
    commandstr=sprintf('%s %s "%s"',path_to_python, path_to_converter, path_to_data);
    system(commandstr);
    clear commandstr fname
    d=dir('tr-*.mat');
    disp(sprintf('%d trace matfiles found\n',length(d))) 
    for c=1:length(d)
        tr=load(d(c).name);
        disp(sprintf('Loaded trace %d',c))
        scnl(c)=scnlobject(tr.station, tr.channel, tr.network, tr.location);
        snum=datenum('2007-08-28T00:00:00','yyyy-mm-ddTHH:MM:SS');
        w(c)=waveform(scnl, tr.sampling_rate, snum, tr.data);
        delete(d(c).name);
    end
end
