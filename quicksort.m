function arr = quicksort(arr, start_sort, end_sort)

    if nargin == 0
        arr = rand(25,1);
        arr = quicksort(arr);
        return
    end

    if nargin < 2
        N = length(arr);
        start_sort = 1;
        end_sort = N;
    end

    if start_sort >= end_sort
        return
    end
    
    [arr,ind] = partition_sort(arr, start_sort, end_sort);
    
    arr = quicksort(arr,start_sort, ind-1);
    arr = quicksort(arr,ind+1,end_sort);
    
end

function [arr,ind] = partition_sort(arr, start_sort, end_sort)
    global h_bar_quicksort;
    if isempty(h_bar_quicksort)
        h_bar_quicksort = bar(arr);
        h_bar_quicksort.FaceColor = 'flat';
        drawnow
    else
        try 
            h_bar_quicksort.YData = arr;
        catch
            h_bar_quicksort = bar(arr);
            h_bar_quicksort.FaceColor = 'flat';
            drawnow
        end            
    end
    pivotIndex = start_sort;
    pivotValue = arr(end_sort);
    for k = start_sort:end_sort-1
        if arr(k) < pivotValue
            arr = swap(arr, k, pivotIndex);
            h_bar_quicksort.YData = arr;
            pivotIndex = pivotIndex + 1;
        end
        h_bar_quicksort.CData = repmat([0,0,1],length(arr),1);
        h_bar_quicksort.CData(start_sort:end_sort,:) = repmat([0,1,0],end_sort-start_sort+1,1);
        h_bar_quicksort.CData(pivotIndex,:) = [1,0,0];
        drawnow
    end
    arr = swap(arr, pivotIndex, end_sort);
    h_bar_quicksort.YData = arr;
    h_bar_quicksort.CData = repmat([0,0,1],length(arr),1);
    h_bar_quicksort.CData(start_sort:end_sort,:) = repmat([0,1,0],end_sort-start_sort+1,1);
    h_bar_quicksort.CData(pivotIndex,:) = [1,0,0];
    drawnow
    ind = pivotIndex;
end



function arr = swap(arr,a,b)

    temp = arr(a);
    arr(a) = arr(b);
    arr(b) = temp;
    
end
        