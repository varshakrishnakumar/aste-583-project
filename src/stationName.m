function name = stationName(id)
    switch id
        case 1
            name = 'Goldstone';
        case 2
            name = 'Canberra';
        case 3
            name = 'Madrid';
        case 4
            name = 'Antarctica';
        otherwise
            name = 'Unknown';
    end
end