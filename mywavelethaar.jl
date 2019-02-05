function wavelethaar(data, stoplevel = 6, ds = 0 ; level = 0)
    if ds == 0 
        tmp = size(data)
        b = tmp[1] > tmp[2] ? tmp[1] : tmp[2]
        ds = 2^Int(ceil(log2(b))) 
    end
    if level >= stoplevel || ds == 1
        return data
    elseif ds == 0
        return nothing
    end
    cdata = copy(data[1:ds, 1:ds])
    hs = ds ÷ 2
    for i in 1:hs
        for j in 1:hs
            tl = data[2i - 1, 2j - 1]
            tr = data[2i, 2j - 1]
            bl = data[2i - 1, 2j]
            br = data[2i, 2j]

            tlrp = tl + tr 
            tlrm = tl - tr
            blrp = bl + br
            blrm = bl - br

            ntl = (tlrp + blrp) / 4
            ntr = (tlrm + blrm) / 4
            nbl = (tlrp - blrp) / 4
            nbr = (tlrm - blrm) / 4

            cdata[2i - 1, 2j - 1] = ntl
            cdata[2i, 2j - 1] = ntr
            cdata[2i - 1, 2j] = nbl
            cdata[2i, 2j] = nbr
        end
    end

    for i in 1:hs
        for j in 1:hs
            data[i, j] = cdata[2i - 1, 2j - 1]
            data[i + hs, j] = cdata[2i, 2j - 1]
            data[i, j + hs] = cdata[2i - 1, 2j]
            data[i + hs, j + hs] = cdata[2i, 2j]
        end
    end
    wavelethaar(data, stoplevel, ds ÷ 2, level = level+1)
end
