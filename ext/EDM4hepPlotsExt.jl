module EDM4hepPlotsExt

    using EDM4hep
    using EDM4hep.Histograms
    using FHist, RecipesBase, Plots, Statistics

    @recipe function f(h::H1D)
        seriestype := :stepbins
        hist = h.hist
        N = nentries(hist)
        M = round(mean(hist); sigdigits= 2)
        S = round(std(hist); sigdigits= 2)
        label --> "Entries = $N\nMean = $M\nStd Dev = $S"
        title := h.title * (h.usym != :nounit ? " [$(h.usym)]" : "")
        xlabel := string(h.usym)
        ylabel := "Entries"
        x:= binedges(hist)
        y:= bincounts(hist)
        ()
    end

    @recipe function f(h::H2D)
        hist = h.hist
        seriestype := :bins2d
         title := h.title
        x := binedges(hist)[1]
        y := binedges(hist)[2]
        z := (surf = bincounts(hist), )
        ()
    end

end