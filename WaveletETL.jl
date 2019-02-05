module WaveletETL

using FileIO
using ETLCDBReader
using Images
using Wavelets
export getfloatimg, mostnear

include("mywavelethaar.jl")

function getfloatimg(img)
    f = (zeros(Float32,64,64))
    f[1:63, 1:64] = img.data
    f
end

getcharetls(cc::Char, etls) = filter(x -> x.charcode == cc, etls)

function heiabs(cc::Char, etls; r = 1:100, level=5, wavefun=WT.haar)
    allcc = getcharetls(cc, etls)
    h = Gray.(zeros(Float32,64,64))
    d = (zeros(Float32,64,64))
    yoko0 = (zeros(Float32, 64))
    for ne in allcc[r]
        d[1:63, 1:64] = ne.data
        d[64, 1:64] = yoko0
        d = wplotim(d, level, wavelet(wavefun))
        h += (d)
    end
    h / length(r)
end

waveabs(img, level=5; wavefun=WT.haar) = Gray.(wplotim(img, level, wavelet(wavefun)))

function kura(charetls::ETLCDBReader.ETL, heiimg, level=5; wavefun=WT.haar)
    a = Float32[]
    d = (zeros(Float32, 64, 64))
    wd = (zeros(Float32, 64, 64))
    for ne in charetls[101:200]
        d[1:63, 1:64] = ne.data
        wd = wplotim(d, level, wavelet(wavefun))
        push!(a, sum(abs, (heiimg - wd)))
    end
    a
end

function kura(charetls::Vector{Matrix{Float32}}, heiimg, level=5; wavefun=WT.haar)
    a = Float32[]
    for ne in charetls
        d = wplotim(ne, level, wavelet(wavefun))
        push!(a, sum(abs ,(heiimg - d)))
    end
    a
end

function kura2(charetls::ETLCDBReader.ETL, heiimg, level=5; wavefun=WT.haar)
    a = Float64[]
    d = (zeros(Float64, 64, 64))
    yoko0 = (zeros(Float64, 64))
    for ne in charetls[101:200]
        d[1:63, 1:64] = ne.data
        d[64, 1:64] = yoko0
        d = wplotim(d, level, wavelet(wavefun))
        res = sum(map((x,y) -> if y == 0 || x == 0 0 else abs(x-y) end , heiimg, d))
        push!(a, res)
    end
    a
end

function mostnear(allhrgn, etls, heis; level=5, wavefun=WT.haar)
    for hrgn in allhrgn
        taretls = WaveletETL.getcharetls(hrgn, etls)
        mhl = Dict(map(x -> x => 0 ,allhrgn))
        tardiffs = Dict{Char, Array{Float32}}()
        for hc in allhrgn
            tardiffs[hc] = kura(taretls, heis[hc], level, wavefun=wavefun)
        end
        for i in 1:100
            m = sort(map(hc -> hc => tardiffs[hc][i], allhrgn), by=x->x.second)[1]
            mhl[m.first] += 1
        end
        println("$(hrgn): mostnear:$(sort(collect(mhl), by=x->x.second, rev=true)[1])")
        println("$(hrgn)と判定される割合: $(mhl[hrgn]/100)")
    end
end

getheis(charlist::Vector{Char}, etls; level=5, wavefun=WT.haar) = Dict(map(x -> x => heiabs(x, etls, level=level, wavefun=wavefun), charlist))

function test()
	etls = ETLCDBReader.getetl9b("/home/pisces/Downloads/ETL9B/")
    heis = getheis(allhrgn, etls)
	mostnear(allhrgn, etls, heis6)
end

function savewaveletimg(char, etls; level=5, wavefun=WT.haar)
    FileIO.save("$(char * level).bmp", waveabs(getfloatimg(getcharetls(char, etls)[1]), level))
end

function main()
    etls = ETLCDBReader.getetl9b("/home/pisces/Downloads/ETL9B")
    d = getfloatimg(etls[39].data)
    return wplotim(d, 5, wavelet(WT.haar))
end

function hoge()
    #for i in 1:64
    #    for j in 1:64
    #        print("$(d[j,i]), ")
    #        #print("\e[48;5;$(0xE8+abs(trunc(Int,23*d[j,i])))m \e[m")
    #        #print("\e[48;5;$(0xE8+abs(trunc(Int,23*d[j,i])))m \e[m")
    #    end
    #    print("\n")
    #end
end

allhrgn = ['あ', 'い', 'う', 'え', 'お', 'か', 'が', 'き', 'ぎ', 'く', 'ぐ', 'け', 'げ', 'こ', 'ご', 'さ', 'ざ', 'し', 'じ', 'す', 'ず', 'せ', 'ぜ', 'そ', 'ぞ', 'た', 'だ', 'ち', 'ぢ', 'つ', 'づ', 'て', 'で', 'と', 'ど', 'な', 'に', 'ぬ', 'ね', 'の', 'は', 'ば', 'ぱ', 'ひ', 'び', 'ぴ', 'ふ', 'ぶ', 'ぷ', 'へ', 'べ', 'ぺ', 'ほ', 'ぼ', 'ぽ', 'ま', 'み', 'む', 'め', 'も', 'や', 'ゆ', 'よ', 'ら', 'り', 'る', 'れ', 'ろ', 'わ', 'を', 'ん']

function checkdiff(cs, etls, l; wavefun=WT.haar)
    a = cs[1]
    b = cs[2]
    alla = getcharetls(a, etls)
    allb = getcharetls(b, etls)
    heia = heiabs(a, etls, level=l, wavefun=wavefun)
    heib = heiabs(b, etls, level=l, wavefun=wavefun)
    aadiffs = kura(alla, heia, l, wavefun=wavefun)
    abdiffs = kura(alla, heib, l, wavefun=wavefun)
    badiffs = kura(allb, heia, l, wavefun=wavefun)
    bbdiffs = kura(allb, heib, l, wavefun=wavefun)
    akaku = sum(map((x,y) -> x < y ? 1 : 0, aadiffs, abdiffs)) / 100
    bkaku = sum(map((x,y) -> x < y ? 1 : 0, bbdiffs, badiffs)) / 100
    println("level:$(l)")
    println("Char:$(a)が$(a)と判定される割合は$(akaku)")
    println("Char:$(b)が$(b)と判定される割合は$(bkaku)")
    (akaku, bkaku)
end

#if length(PROGRAM_FILE)!=0 && contains(@__FILE__, PROGRAM_FILE)
#main()
#end

end

# めも
#
# 前期の活動
# 予備実験
# 展望 (今後の計画)
# 引用文献
