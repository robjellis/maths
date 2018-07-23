function data = tapread( filename )

fid = fopen( filename );

if fid == -1
    ['problem opening text file: ' filename]
else
    % read the data. the timestamps are in c{1}
    c = textscan( fid, '%f' );

    data = c{1};
    
    fclose( fid );
end