package picard.vcf.filter;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.HashMap;
import java.util.Map;

/**
 * Filters out a record if all variant samples have depth lower than the given value(s).
 */
public class DepthFilter implements GenotypeFilter {
    private final int minHetSnp, minHomSnp, minHetIndel, minHomIndel;

    public DepthFilter(final int minHetSnp, final int minHomSnp, final int minHetIndel, final int minHomIndel) {
        this.minHetSnp = minHetSnp;
        this.minHomSnp = minHomSnp;
        this.minHetIndel = minHetIndel;
        this.minHomIndel = minHomIndel;
    }

    @Override
    public String filter(final VariantContext ctx, final Genotype gt) {
        final int minHet, minHom;
        if (ctx.isSNP()) {
            minHet = minHetSnp;
            minHom = minHomSnp;
        }
        else {
            minHet = minHetIndel;
            minHom = minHomIndel;
        }

        if (!gt.isHomRef() && gt.getDP() < (gt.isHet() ? minHet : minHom)) {
            return "LowDP";
        }
        else {
            return null;
        }
   }
}
