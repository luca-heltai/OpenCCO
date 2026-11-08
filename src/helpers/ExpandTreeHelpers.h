/**
 *  OpenCCO implementation 
 *  Copyright (C) 2023 B. Kerautret;  Phuc Ngo, N. Passat H. Talbot and C. Jaquet
 *
 *  This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 **/

#pragma once

#if defined(CONSTRUCTIONHELPERS_RECURSES)
#error Recursive header files inclusion detected in ConstructionHelpers.h
#else // defined(CONSTRUCTIONHELPERS_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CONSTRUCTIONHELPERS_RECURSES

#if !defined CONSTRUCTIONHELPERS_h
/** Prevents repeated inclusion of headers. */
#define CONSTRUCTIONHELPERS_h

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/io/readers/GenericReader.h"
#include "CoronaryArteryTree.h"



namespace ExpandTreeHelpers {


template<typename DomCtr, int TDim>
inline
void
initFirtElemTree(CoronaryArteryTree< DomCtr, TDim > &aTree,
                 bool verbose = false){
    for(unsigned int i = 0; i < TDim; i++){
        aTree.myTreeCenter[i] = aTree.myDomainController().myCenter[i];
    }
    
    auto p = aTree.myDomainController().firtCandidatePoint();
    aTree.myVectSegments[0].myCoordinate = p;
}


// Insert one explicit site through the same native optimizer and acceptance
// gate as sampled growth. Proposal spacing is a sampling policy; lumen
// clearance and all native/caller geometry checks still apply here.
template<typename DomCtr, int TDim, typename AcceptCandidate>
inline bool
tryAddPoint(CoronaryArteryTree<DomCtr,TDim> &tree,
            const typename CoronaryArteryTree<DomCtr,TDim>::TPointD &point,
            bool verbose, AcceptCandidate acceptCandidate)
{
    if (tree.myKTerm >= tree.my_NTerm) return false;
    for (unsigned int j=0; j<TDim; ++j)
        if (!std::isfinite(point[j])) return false;
    for (unsigned int i=1; i<tree.myVectSegments.size(); ++i)
        if (!(tree.getProjDistance(i,point) > tree.myVectSegments[i].myRadius)) return false;
    bool found=false;
    auto best=tree;
    double bestVolume=0.;
    for (const auto index : tree.getN_NearestSegments(point,tree.myNumNeighbor)) {
        const auto center=tree.findBarycenter(point,index);
        if (tree.isIntersectingTree(point,center,tree.myVectSegments[index].myRadius,index)) {
            if (tree.myCandidateDiagnostics) ++(*tree.myCandidateDiagnostics)["prefilter_collision"];
            continue;
        }
        auto candidate=tree;
        if (!candidate.isAddable(point,index,100,.01,candidate.myNumNeighbor,verbose)) continue;
        candidate.updateLengthFactor();
        candidate.updateResistanceFromRoot();
        candidate.updateRootRadius();
        if (!acceptCandidate(candidate)) continue;
        const double volume=candidate.computeTotalVolume(1);
        if (!found || volume<bestVolume) {
            bestVolume=volume;
            best=candidate;
        }
        found=true;
    }
    if (found) tree=best;
    return found;
}

template<typename DomCtr, int TDim, typename AcceptCandidate>
inline bool
expandTree(CoronaryArteryTree<DomCtr,TDim> &tree,
           bool verbose, unsigned int nbMaxSearch,
           unsigned int nbTryCandidate, AcceptCandidate acceptCandidate)
{
    // The caller owns the RNG seed. Never restart it from wall-clock time here.
    while (tree.myKTerm<tree.my_NTerm) {
        DGtal::trace.progressBar(tree.myKTerm,tree.my_NTerm);
        bool found=false;
        for (unsigned int search=0; search<nbMaxSearch && !found; ++search) {
            const auto sample=tree.tryGenerateNewLocation(nbTryCandidate);
            if (!sample.second) {
                if (tree.myCandidateDiagnostics) ++(*tree.myCandidateDiagnostics)["sampling_budget"];
                continue;
            }
            found=tryAddPoint(tree,sample.first,verbose,acceptCandidate);
        }
        if (!found) break;
    }
    const bool complete=tree.myKTerm==tree.my_NTerm;
    if (!complete)
        DGtal::trace.warning()<<"CCO budget exhausted: partial tree ("
                              <<tree.myKTerm<<" of "<<tree.my_NTerm
                              <<" requested terminals)"<<std::endl;
    return complete;
}

// Legacy callers retain the same candidate selection and final normalization.
template<typename DomCtr, int TDim>
inline bool expandTree(CoronaryArteryTree<DomCtr, TDim> &tree,
                       bool verbose = false, unsigned int searches = 100,
                       unsigned int samples = 1000) {
    return expandTree(tree, verbose, searches, samples, [](const auto &) { return true; });
}

template<typename TImage >
std::vector<std::vector<typename TImage::Domain::Point > >
getImageContours(const TImage &image,
                 unsigned int threshold=128){
    std::vector<std::vector<typename TImage::Domain::Point > > v;
    DGtal::trace.error() << "Use CCO is only implemented in 2D and 3D and use ImageContainerBySTLVector with unsigned char"
    << "to export domain."
    << "You use such an image : " << image <<  std::endl;
    throw 1;
    return v;
}


std::vector<std::vector<DGtal::Z2i::Point > >
getImageContours(const DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> &image,
                 unsigned int threshold){
    typedef  DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> TImage;
    typedef DGtal::functors::IntervalThresholder<typename TImage::Value> Binarizer;
    DGtal::Z2i::KSpace ks;
    if(! ks.init( image.domain().lowerBound(),
                 image.domain().upperBound(), true )){
        DGtal::trace.error() << "Problem in KSpace initialisation"<< std::endl;
    }
    
    Binarizer b(threshold, 255);
    DGtal::functors::PointFunctorPredicate<TImage,Binarizer> predicate(image, b);
    DGtal::trace.info() << "DGtal contour extraction from thresholds ["<<  threshold << "," << 255 << "]" ;
    DGtal::SurfelAdjacency<2> sAdj( true );
    std::vector< std::vector< DGtal::Z2i::Point >  >  vectContoursBdryPointels;
    DGtal::Surfaces<DGtal::Z2i::KSpace>::extractAllPointContours4C( vectContoursBdryPointels, ks, predicate, sAdj );
    return vectContoursBdryPointels;
}




/**
 Template Specialisation in 3D  to export surfel image  border of the restricted image domain.
 * Todo @BK
 */
std::vector<std::vector<DGtal::Z3i::Point > >
getImageContours(const DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, unsigned char> &image,
                 unsigned int threshold){
    std::vector<std::vector<DGtal::Z3i::Point > > v;
    
    return v;
}



}

#endif // !defined CONSTRUCTIONHELPERS_h

#undef CONSTRUCTIONHELPERS_RECURSES
#endif // else defined(CONSTRUCTIONHELPERS_RECURSES)




