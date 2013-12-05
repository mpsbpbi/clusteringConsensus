#!/usr/bin/env python
import pdb
import os
import sys
import random
import itertools as it
import re
import subprocess
import logging
import time
import h5py

from optparse import OptionParser

from pbpy.align.Exonerate import *
from pbpy.align.Lastz import *
from pbpy.align.Blasr import *
from pbpy.align.YassAlign import *
from pbpy.align.Bwasw import *
from pbpy.align.ZScore import *
from pbpy.model.Range import Range, Ranges
from pbpy.util.Process import backticks
from pbpy.util.Timing import *
from pbpy.io.FastaIO import SimpleFastaReader
from pbpy.io.GffIO import GffReader
from pbpy.io.ReferenceEntry import ReferenceEntry
from pbpy.io.VersionUtils import getSmrtAnalysisVersionInfo

from pbpy.io import cmph5
from pbpy.io.cmph5.CmpH5 import SFCmpH5
from pbpy.io.cmph5.CmpH5 import PBCmpH5
from pbpy.io.cmph5.CmpH5Utils import PBCmpH5MetaData, repositoryPathToURI
from pbpy.io.MovieHDF5IO import getSoftwareVersion


# defines a large reference for the purpose of more aggressive splitting of the
# exonerate jobs
LARGE_REFERENCE_LENGTH = 100000
#LARGE_REFERENCE_LENGTH = 200000
EMPTY_READ_GROUP = 'NULL_GROUP'

log = logging.getLogger(__name__)


class CompareSequences(object):
    """Compares sequences to each other using an algorithm
    selected from a selection of supported command-line
    alignment algorithms.
    Output is in cmp.h5.

    Miscellaneous options for parallel processing,
    calculating read significance and filtering alternative hits.
    """
    VERSION = 1.1

    def __init__(self, argv):
        self._tmpRepos = None
        self.__parseOptions(argv)
        self._refSeqMap = {}
        if self.options.debug:
            self.startMemory = memory()
            self.startResident = resident()


    def _logMemory( self, msg ):
        memUsed = memory( self.startMemory )
        resUsed = resident( self.startResident )
        print >>sys.stderr, "%s, virtual memory used (MB): %.1f" % (msg,memUsed/1024.0/1024.0)
        print >>sys.stderr, "%s, resident memory used (MB): %.1f" % (msg,resUsed/1024.0/1024.0)

    def __parseOptions( self, argv ):
        usage = 'Usage: %prog [--help] [options] seq1_query.fsta seq2_target.fsta > alignment.out'
        parser = OptionParser( usage=usage, description=CompareSequences.__doc__ )
        # pass exonerate options by prepending the option with '-x '
        # i.e. to specify --seedrepeat 3 use -x --seedrepeat 3
        parser.add_option( '--dryrun', action="store_true", help="Skip alignment")
        parser.add_option( '--debug', action="store_true",
            help="Turn on debugging output" )
        parser.add_option( '--info', action="store_true",
            help="Turn on basic output" )
        parser.add_option( '--algorithm',
            help="Alignment algorithm [exonerate|yass|blasr|lastz]" )
        parser.add_option( "-x", nargs=2, dest='algorithmOpts',
            action="append", help='Pass through options for algorithm (e.g. exonerate)' )
        parser.add_option( "--local", action="store_true",
            help="Use local alignment instead of global alignment" )
        parser.add_option( "--dirty", action="store_true",
            help="Use settings that are good for aligning dirty reads to a reference" )
        parser.add_option( "--minAccuracy", type="float",
            help="Min accuracy to output a hit" )
        parser.add_option( "--minLength", type="float",
            help="Min length to output a hit" )
        parser.add_option( "--minZ", type="float",
            help="Min Z score to output a hit")
        parser.add_option( "--uniqueIds", action="store_true",
            help="Modify query ids so that each hit has a unique id" )
        parser.add_option( "--nproc", type="int",
            help="Number of processors to use for alignment computation" )
        parser.add_option( "--delta", action="store_true",
            help="Output hits in the delta format from MUMmer" )
        parser.add_option( '--noiseData',
            help="noise triplet, .xy or compare XML file containing noise data for estimating a z-score" )
        parser.add_option( "--trimWindow", type="int",
            help="Window size for trimming ends [not implemented]" )
        parser.add_option( "--trimErrors", type="int",
            help="Trim ends until there are less than trimErrors in the window [not implemented]" )
        parser.add_option( "--multiple",
            help="Specify a policy for how to treat multiple hits. One of [random,all,deltaz,bestscore,leftpositive]\n bestscore (DEFAULT) returns the best alignment score hit.  All returns all hits, limited by the alignment routine, and deltaz limits by z value." )
        parser.add_option( "--h5fn",
            help="Write to specified file in HDF5 format, set showAlignmet to be true " )
        parser.add_option( "--h5mode",
            help=" 'w' for creating a new hdf5 file, 'a' for appending " )
        parser.add_option( "--h5pbi", action="store_true",
            help="Create a PBCmpH5 file. Only works with .fofn input and the --h5fn option." )
        parser.add_option( "--refSeqName",
            help=" string for the name of the target sequence " )
        parser.add_option( "--useTemp", dest='manualTempFileName',
            help=" Specify a temporary output name, rather than using an auto-generated one in/tmp")
        parser.add_option( "--tmpDir", dest='tmpDir', help="temporary output dir for blasr", default="/scratch")
        parser.add_option( "--pls2fasta", action="store_true",
            help=" Convert pls files into fasta before alignment.")
        parser.add_option( "--noSplitSubreads", action="store_true",
            help=" Do not split reads into subreads even if subread regions are available.")
        parser.add_option( "--advanceExactMatches",
            help=" Speed up alignments in BLASR by not finding matches of subsequences that have already been matched.  The method skips the last L-n bases of a match, where the match is of length L, and n is specified as the value for advanceExactMatches.")
        parser.add_option( "--advanceHalf", action="store_true",
            help=" Hack for speeding up BLASR. This is deprecated as of v1.3.1")
        parser.add_option( "--forPicard", action="store_true",
            help=" Modify output so that it is compatible with Picard parsing: interleaving ins/del/ins are reordered to be grouped by ins and del.")
        parser.add_option( "--regionTable",
            help=" Specify a regions table for filtering reads.")
        parser.add_option("--useQuality", action="store_true",
            help=" Force BLASR to use quality values when computing a local alignment.")
        parser.add_option("--ignoreRegions", action="store_true", \
            help=" Ignore all information in a regions table, even if one exists")
        parser.add_option("--ignoreHQRegions", action="store_true", \
            help=" Do not use high-quality region information, even if it exists (no masking is done).")
        parser.add_option( "--noXML", action="store_true", \
            help="turning off XML generation. Deprecated. XML Output is no longer supported. " )
        parser.add_option( "--singleReadGroup", action="store_true", \
            help="When writing cmpH5 store reads in a single default read group." )
        parser.add_option( "--subsample", type = "float", \
            help="Makes BLASR subsample and align reads at a desired fraction of input." )
        parser.add_option( "--randomizeDeletions", action="store_true",
            help="Post-process alignment using a gap randomizer for deletions in homopolymers" )
        parser.add_option( "--filterAdapterOnly", action="store_true",
            help="If specified, do not report adapter-only hits using annotations associated with the reference entry" )
        parser.add_option( "--respectFastaGivenSubreadLocation", action="store_true",
            help="If specified, add the given subread start to the hit coordinates of the query_id in the cmp.h5 file giving global read location" )
        parser.add_option( "--useGuidedAlign", action="store_true",
            help="Pass the useGuidedAlign option to BLASR" )
        parser.add_option( "--extend", action="store_true",
            help="Pass the extend option to BLASR.  Extend alignments past sparse dynamic programming hits" )
        parser.add_option( "--divideByAdapter", action="store_true",
            help="Pass the divideByAdapter option to BLASR. Divide reads into subreads using the adapter sequence." )
        parser.add_option( "--useCcs",
            help="[fullpass|allpass|denovo].  Map the ccsSequence to the genome first, then align subreads to the interval that the CCS read mapped to.  fullpass only aligns subreads that span the length of the template.  Specifying allpass maps all subreads.")
        parser.add_option( "--placeRepeatsRandomly", action="store_true",
            help="When there are multiple positions to map a read with equal alignment scores, place a read randomly at one of them.  The default is to place the read at the first.")
        parser.add_option("--concordant", action="store_true",
            help="Map all subreads of a zmw (hole) to where the longest full pass subread of the zmw aligned to.")
        parser.add_option( "--seed", type='int',
            help="Initialize the random number generator with a none-zero integer. 0 means current system time is used.")
        parser.add_option( "--forwardOnly", action="store_true",
            help="Map only to the forward strand.")
        parser.add_option( "--profile", action="store_true",
            help="Use the python cProfile module to profile the running of this script" )
        parser.add_option( "--readType",
            help="Specify the 'ReadType' attribute in the cmp.h5 output")
        parser.add_option( "--nucmer",
            help="Make blasr emulate nucmer (for mapping very long contigs)")
        parser.add_option( "--scoreType",
                           help = "Inform BLASR which method to use for scoring alignments: 0 = edit distance, 1 = quality value sum")
        parser.add_option( "--circularReference", action="store_true",
            help="Map to circular contig. Replicate reference sequence into a concatamer. Collapse hits' coordinates and break them up across repeat boundaries. (BLASR only)")


        parser.set_defaults(dryrun=False,
                            debug=False,
                            info=False,
                            algorithm="exonerate",
                            local=False,
                            dirty=False,
                            nproc=1,
                            delta=False,
                            showAlignment=False,
                            noiseData=None,
                            trimWindow=0,
                            trimErrors=0,
                            minAccuracy=0.0,
                            uniqueIds=False,
                            minLength=0,
                            minZ=-1000000,
                            multiple="bestscore",
                            h5mode='w',
                            h5pbi=False,
                            noXML=False,
                            singleReadGroup=False,
                            refSeqName="ref",
                            regionTable="",
                            manualTempFileName="",
                            #ignoreQualities=False,
                            noSplitSubreads=False,
                            advanceExactMatches=0,
                            advanceHalf=False,
                            forPicard=False,
                            pls2fasta=False,
                            subsample=1.0,
                            randomizeDeletions=False,
                            filterAdapterOnly=False,
                            respectFastaGivenSubreadLocation=False,
                            useGuidedAlign=False,
                            useCcs="",
                            profile=False,
                            readType="standard",
                            circularReference=False,
                            seed=0)

        self.options, self.args = parser.parse_args( argv )

        if not os.path.exists(self.options.tmpDir):
            os.makedirs(self.options.tmpDir)


        #set up the logger now
        if self.options.debug:
            level = logging.DEBUG
        else:
            level = logging.INFO
        self._setupLog(level)

        #initial the random generator if --seed is set
        if self.options.seed != 0:
            random.seed(self.options.seed)
            log.debug("Set seed to {0}".format(self.options.seed))

#        for i in range(0,len(self.options)):
#            print "option " + i + " " + self.options[i]
        if self.options.h5fn != "":
            self._outputHDF5 = True
            self.options.showAlignment = True
        else:
            self._outputHDF5 = False

        if len(self.args)!=3:
            parser.error( 'Expected 2 arguments' )

        # escape spaces since these file names are passed off to processes
        self.sequence1 = self.args[1].replace( ' ', r'\ ' )

        #guarantees that the reference is from refRepos
        self.sequence2 = self._getReferenceDir( self.args[2].replace( ' ', r'\ ' ) )
        self.reference = ReferenceEntry( self.sequence2 )

        if self.options.algorithmOpts:
            self.algorithmParameters = dict( self.options.algorithmOpts )
        else:
            self.algorithmParameters = {}

        if self.options.pls2fasta:
            #
            # Convert pls files into temporary fasta files, and create a new fofn if the current one is a fofn.
            #

            self.tmpFastaFileName = tempfile.mkstemp(dir=self.options.tmpDir, suffix=".fasta")[1]
            opts = []
            if (self.options.regionTable != ""):
                opts += [ "-regionTable", self.options.regionTable, "-trimByRegion" ]
            if (self.options.noSplitSubreads == True):
                opts += [ "-noSplitSubreads" ]

            pls2fastaCmdLine = "pls2fasta " + self.sequence1 + " " + self.tmpFastaFileName + " " + " ".join( opts )

            pls2fastaProc = subprocess.Popen( pls2fastaCmdLine, shell=True )
            retCode = pls2fastaProc.wait()
            self.sequence1 = self.tmpFastaFileName

        if self.options.h5pbi and self.options.noiseData is None:
            parser.error( 'Use of --h5pbi requires --noiseData as '
                          'well ( e.g. --noiseData="-77.27,0.08654,0.00121" ).' )

        if self.options.algorithm == "exonerate":
            # exonerate specific options
            self.__setDefaultExonerateParameters()
            if self.options.local:
                self.algorithmParameters['-m'] = 'affine:local'
            if self.options.dirty:
                self.algorithmParameters['-m'] = 'affine:local'
                self.algorithmParameters['--gapopen']   = '-7'
                self.algorithmParameters['--gapextend'] = '-9'
                self.algorithmParameters['--percent']   = '5'
        elif self.options.algorithm == "blasr":
            # blasr specific options
            # the -m option is the BLASR output format
            # 5 is required for parsing by BlasrService
            self.algorithmParameters['-m'] = '5'
            if self.options.regionTable != "":
                self.algorithmParameters.setdefault( '-regionTable', self.options.regionTable)

            self.algorithmParameters.setdefault( '-maxExpand', '1' )
            self.algorithmParameters.setdefault( '-maxScore', '-100' )
            if self.options.noSplitSubreads == True:
                self.algorithmParameters['-noSplitSubreads'] = True
            if self.options.nucmer == True:
                self.algorithmParameters['-nucmer'] = True
            if self.options.advanceHalf == True:
                self.algorithmParameters['-advanceHalf'] = True
            if self.options.forPicard == True:
                self.algorithmParameters['-forPicard'] = True
            if self.options.advanceExactMatches != 0:
                self.algorithmParameters['-advanceExactMatches'] = self.options.advanceExactMatches
            if self.options.useQuality == True:
                self.algorithmParameters['-useQuality'] = ' '
            if self.options.ignoreRegions == True:
                self.algorithmParameters['-ignoreRegions'] = True
            if self.options.ignoreHQRegions == True:
                self.algorithmParameters['-ignoreHQRegions'] = True
            if self.options.useGuidedAlign == True:
                self.algorithmParameters['-nouseDetailedSDP'] = True
                self.algorithmParameters['-useGuidedAlign'] = True
            if self.options.subsample < 1.0:
                self.algorithmParameters['-subsample'] = self.options.subsample
            if self.options.useCcs == 'fullpass':
                self.algorithmParameters['-useccs'] = True
            if self.options.placeRepeatsRandomly == True:
                self.algorithmParameters['-placeRepeatsRandomly'] = True
            if self.options.concordant == True:
                self.algorithmParameters['-concordant'] = True
            if self.options.seed != 0:
               self.algorithmParameters['-randomSeed'] = self.options.seed;
               print "self.algorithmParameters[-randomSeed] = {0}\n".format(self.options.seed)
            if self.options.forwardOnly == True:
                self.algorithmParameters['-forwardOnly'] = True
            if self.options.useCcs == 'allpass':
                self.algorithmParameters['-useccsall'] = True
            if self.options.useCcs == 'denovo':
                self.algorithmParameters['-useccsdenovo'] = True
            if self.options.scoreType != None:
                self.algorithmParameters['-scoreType'] = self.options.scoreType
            if self.options.extend == True:
                self.algorithmParameters['-extend'] = True
            if self.options.divideByAdapter == True:
                self.algorithmParameters['-divideByAdapter'] = True

        elif self.options.algorithm == "lastz":
            # lastz works pretty well by default
            self.algorithmParameters['--format']='=lav+txt'
        else:
            # Set any other algorithm specific options here.
            pass
        # resolve incompatible options
        self._makeSane()

    def _makeSane(self):
        '''
        Resolve incompatible options and return error message.
        '''
        if (self.options.placeRepeatsRandomly == True and self.options.multiple != "bestscore"):
            msg = "PlaceRepeatsRandomly can not be used with multiple hit policy [{0}]"\
                  .format(self.options.multiple)
            log.error(msg)
            raise SystemExit(msg)


    def _getReferenceDir(self, refPath):
        '''
        A ReferenceEntry must be instantiated from a reference directory.
        If refPath is a reference dir, return it.
        If a fasta file is passed in, create a reference using referenceUploader and return
        the resulting reference dir.
        '''
        if os.path.isdir(refPath):
            return refPath

        #Allow user to easily identify temp repos dir
        suff = os.path.basename(refPath)
        try:
            suff = suff.split()[0][:5]
        except:
            pass

        self._tmpRepos = tempfile.mkdtemp(suffix=suff, prefix="repos", dir=self.options.tmpDir);

        #hard-code the name, so we're guaranteed of legal chars
        name = "reference"
        res = backticks( 'referenceUploader -c -n%s -f%s -p%s' % (name, refPath, self._tmpRepos) )
        if res[1] != 0:
            raise SystemExit, 'Failed to create temp reference repository: %s' % res[2]

        refDir = os.path.join( self._tmpRepos, name )
        log.debug("Invoked referenceUploader to create temp reference: %s" % refDir )
        return refDir

    def __setDefaultExonerateParameters(self):
        if '-m' not in self.algorithmParameters:
            self.algorithmParameters[ '-m' ] = 'affine:global'
        if '--exhaustive' not in self.algorithmParameters:
            self.algorithmParameters[ '--exhaustive' ] = 'TRUE'
        # if the user specified delta-z filtering of alternate hits
        # then we override the defaults and enforce allowing alternate
        # hits (seems reasonable)
        if 'deltaz' in self.options.multiple or \
            'all' in self.options.multiple or \
            'leftpositive' in self.options.multiple or \
            'leftpositivelong' in self.options.multiple:
            self.algorithmParameters[ '--bestn' ] = '10'
            self.algorithmParameters[ '--subopt' ] = 'TRUE'
        else:
            if '--bestn' not in self.algorithmParameters:
                self.algorithmParameters[ '--bestn' ] = '1'
            if '--subopt' not in self.algorithmParameters:
                self.algorithmParameters[ '--subopt' ] = 'FALSE'
        if '--percent' not in self.algorithmParameters:
            self.algorithmParameters[ '--percent' ] = '20'
        if '--fsmmemory' not in self.algorithmParameters:
            self.algorithmParameters[ '--fsmmemory' ] = '128'

    def _summarizeHit(self, hit, full=False):
        if self.options.algorithm == "blasr" and \
            hit.query_strand=='-' and hit.target_strand=='+':
            hit.target_strand, hit.query_strand = '-', '+'
            hit.target_end, hit.target_start = hit.target_start, hit.target_end

        al = abs( hit.query_end - hit.query_start )
        nCorrect = al - hit.nIns - hit.nDel - hit.nMismatch
        print '<hit name="%s" unalignedLength="%d" ' \
            'start="%d" end="%d" strand="%s" target_id="%s" targetStart="%d" targetEnd="%d" targetStrand="%s">' % \
            ( hit.query_id, hit.query_length, hit.query_start, hit.query_end, \
              hit.query_strand, hit.target_id, hit.target_start, hit.target_end, hit.target_strand )
        if self.zcalculator:
            z = calculateZ( self.zcalculator, hit )
            print '<zScore value="%.3f"/>' % z
        print '<score value="%d" />' % (hit.score)
        print '<nInsert value="%d" percent="%.2f" />' % ( hit.nIns, pct(hit.nIns,al) )
        print '<nDelete value="%d" percent="%.2f" />' % ( hit.nDel, pct(hit.nDel,al) )
        print '<nMismatch value="%d" percent="%.2f" />' % ( hit.nMismatch, pct(hit.nMismatch, al) )
        print '<nCorrect value="%d" percent="%.2f" />' % ( nCorrect, pct(nCorrect,al) )
        if full:
            print '<alignment><query>'
            print hit.alignedQuery
            print '</query><target>'
            print hit.alignedTarget
            print '</target></alignment>'
        print '</hit>'


    def _saveAlignmentInHDF5(self, hits, commandLine, mode='w'):
        """Newer method for saving hits to cmp.h5 using OO interface and
        v1.0 semantics"""

        #hits = list(hits)
        if self.options.algorithm=='blasr':
            infile = open( self.algorithmParameters['-titleTable'], 'r' )
            get_prefix = lambda x : x.split('|')[0] if x.startswith('ref0') and '|' in x else x
            idMap = dict( [ ( str(i), get_prefix(name) ) for i, name in enumerate(it.imap(lambda l:l.strip(), infile )) ] )
            infile.close()
            idToInternalId = lambda x: idMap[x]
        else:
            idToInternalId = lambda x: x.split("|")[0] if ( "|" in x and x.startswith("ref") ) else x

        if self.options.noiseData:
            self.zcalculator = ZScoreCalculator()
            self.zcalculator.loadNoiseData( self.options.noiseData )
            self.multipleHitPolicy.setZCalculator(self.zcalculator)
        else:
            self.zcalculator = None

        filename = self.options.h5fn
        cmpType = PBCmpH5 if self.options.h5pbi else SFCmpH5
        readType = "CCS" if self.options.useCcs == "denovo" else self.options.readType
        h5f = cmph5.factory.create( filename, mode, cmpType=cmpType, readType=readType )

        if self.options.debug:
            self._logMemory( 'Before output' )

        for contig in self.reference.contigs:
            log.debug("Adding reference %s to cmp.h5" % str(contig.header) )
            h5f.addReference(   contig.header, md5=contig.digest,
                                length=contig.length )

        primaryVersion = None
        if self.sequence1.endswith(".fofn"):
            for line in open( self.sequence1, 'r' ):
                if os.path.exists( line.strip() ):
                    pv = getSoftwareVersion( line.strip() )[1]
                    if primaryVersion == None or pv == primaryVersion:
                        primaryVersion = pv
                    else:
                        log.warning("Combining data from multiple versions of Primary Analysis "
                                "into a single cmp.h5: '%s' and '%s'. Putting just '%s' in the cmp.h5."
                                % (primaryVersion, pv, primaryVersion))
        if primaryVersion != None:
            h5f.attrs["PrimaryVersion"] = primaryVersion

        pbi = None
        if self.options.h5pbi:
            pbi = PBCmpH5MetaData( self.sequence1 )
            h5f.attrs["ReportsFolder"]   = pbi.reportsFolder
            h5f.attrs["PrimaryPipeline"] = pbi.primaryPipeline
            h5f.attrs["Repository"] = repositoryPathToURI( self.sequence2 )
        for line in open( self.sequence1, 'r' ):
            if '.pls.h5' in line or '.bas.h5' in line:
                movieFile = h5py.File(line.strip(),'r')
                mns = movieFile["/ScanData/RunInfo"].attrs["MovieName"]
                ## different Astro/SF file format.
                if isinstance(mns, str):
                    mn = mns
                else:
                    mn = mns[0]
                movieFile.close()

                if self.options.h5pbi:
                    h5f.addMovie( mn, pbi.exps[mn], pbi.runs[mn] )
                else:
                    h5f.addMovie( mn )

        def query2movie( queryId ):
            return queryId.split("/")[0]

        def processHit( hit ):
            """Adjusts things like strand conventions, z-score, etc. Also syncs with refSeqMap"""
            if self.options.trimWindow:
                hit.trimErrors( self.options.trimErrors, self.options.trimWindow )
            if self.options.randomizeDeletions:
                hit.alignedQuery, hit.alignedTarget = \
                    randomize_deletions( hit.alignedQuery, hit.alignedTarget )
            if self.zcalculator:
                z = calculateZ( self.zcalculator, hit )
            else:
                z = None

            # handle RC differences between blasr and exonerate
            if self.options.algorithm == "blasr" and \
                hit.query_strand=='-' and hit.target_strand=='+':
                    hit.target_strand, hit.query_strand = '-', '+'
                    hit.target_end, hit.target_start = hit.target_start, hit.target_end

            if self.options.algorithm == "exonerate" and not self.options.noSplitSubreads:
                if (hit.query_id.find("/") != -1):
                    subStart = int(hit.query_id.split("/")[-1].split("_")[0])
                else:
                    subStart = 0
                hit.query_start += subStart
                hit.query_end   += subStart

            if self.options.algorithm == "bwasw":
                z = 0.0

            # handle fasta specification of subread locations
            if self.options.respectFastaGivenSubreadLocation:
                if (hit.query_id.find("/") != -1):
                    subStart = int(hit.query_id.split("/")[-2].split("_")[0])
                else:
                    subStart = 0
                hit.query_start += subStart
                hit.query_end   += subStart

            hit.zScore = z

            # Map Astro hit names to their SF equivalent
            if re.search("x\d+_y\d+", hit.query_id):
                hit.query_id = cmph5.astro_id_to_springfield_id( hit.query_id )

            if hit.target_id not in self._refSeqMap:
                internalId = idToInternalId( hit.target_id )
                contig = self.reference.getContig( internalId )
                if contig:
                    #self._refSeqMap[ hit.target_id ] = ( contig.displayName, contig.id )
                    self._refSeqMap[ hit.target_id ] = ( contig.header, contig.id )
                else:
                    h5f.addReference( internalId, md5=None,
                                      length=-1 )
                    self._refSeqMap[ hit.target_id ] = ( internalId, internalId )
            return hit

        def byX( hitList, fn ):
            """Given a streaming hit list, returns a streaming ( fn(hit), correspondingHits ) list."""
            currentValue = None
            currentHits = set( )
            for hit in hitList:
                if currentValue == None:
                    currentValue = fn( hit )
                if fn( hit ) == currentValue:
                    currentHits.add( hit )
                else:
                    yield currentValue, currentHits
                    currentValue = fn( hit )
                    currentHits = set([ hit ])
            yield currentValue, currentHits

        def byQuery( hitList ):
            """Given a streaming hit list, returns a streaming ( queryId, correspondingHits ) list."""
            return byX( hitList, lambda h: h.query_id )

        def byMovie( hitList ):
            """Given a streaming hit list, returns a streaming ( movieId, correspondingHits ) list."""
            return byX( hitList, lambda h: query2movie(h.query_id) )

        def validHits( hitList ):
            """Prunes hits according to the multiple hit policy. Processes resulting hits before yielding them.
            Yields only hits that satisfy the hit criteria."""
            #nQ = 0
            for queryId, hits in byQuery( hitList ):
                #if nQ % 1000 == 0:
                #    subprocess.check_call("echo '%s %d ' >> /home/NANOFLUIDICS/dwebster/mem.txt" % ( queryId.split("/")[0], nQ ), shell=True)
                #    subprocess.check_call("ps aux | grep compareSequences | grep -v grep | cut -d ' ' -f '6,7,8' >> /home/NANOFLUIDICS/dwebster/mem.txt", shell=True)
                #nQ += 1
                for hit in self.multipleHitPolicy.handleMultipleHits( list(hits) ):
                    hit = processHit( hit )
                    if self._satisfiesHitCriteria( hit ):
                        yield hit

        self.reference.cacheLookups( ) # Enable fast lookup for reference information
        for movie, hits in byMovie( validHits( hits ) ):
            hitsByRef = { }
            for hit in hits:
                hitsByRef.setdefault( self._refSeqMap[ hit.target_id ], set( ) ).add( hit )
            for ( contigName, contigId ), refHits in hitsByRef.iteritems():
                readgroup = ( "/"+contigId, contigName,  movie )
                self._writeHDF5Hits( refHits, readgroup, h5f, pbi )

        v1, v2 = getSmrtAnalysisVersionInfo()
        h5f.log( "compareSequences.py",
                 v2,
                 str(datetime.datetime.utcnow().isoformat()),
                 " ".join(sys.argv),
                 "Initial Creation" )
        h5f.log( "blasr",
                 v2,
                 str(datetime.datetime.utcnow().isoformat()),
                 commandLine,
                 "Initial Creation" )
        h5f.close()

        if self.options.debug:
            self._logMemory( 'After output' )

    def _writeHDF5Hits( self, hits, readgroup, h5f, pbi ):
        """Given a mapping of ( refGroupPath, refName, movieName ) => [ hits ], writes the
        hits efficiently to the cmp.h5."""
        rfgp, rf, mn = readgroup
        rg = mn
        log.debug("Writing %d hits for key (%s,%s,%s,%s)" % \
                            ( len(hits), rg, mn, rf, rfgp ) )

        if self.options.h5pbi:
            h5f.writeAlignments( hits, rf, rg, mn, pbi.exps[mn], pbi.runs[mn], refGroupPath=rfgp )
        else:
            h5f.writeAlignments( hits, rf, rg, mn, refGroupPath=rfgp )

    def _satisfiesHitCriteria(self, hit):
        """Checks that the hit is above requested thresholds for accuracy
        and hit length.
        Also checks that hit doesn't match adapter only if requested.
        """
        al = hit.nMatch + hit.nIns + hit.nDel + hit.nMismatch
        nCorrect = al - hit.nIns - hit.nDel - hit.nMismatch

        accuracy = nCorrect / float(al)
        z = calculateZ( self.zcalculator, hit )
        passFilter = (accuracy > float(self.options.minAccuracy)) and \
            (al > float(self.options.minLength)) and \
            (z > float(self.options.minZ))
        if not passFilter:
            return False
        if self.options.filterAdapterOnly:
            adapterHit = self._adapters.checkAdapterOnly( hit )
            if adapterHit:
                return False
        return True

    def _printDeltaOutput(self, hits):
        """
        The delta format is documented at http://mummer.sourceforge.net/manual/
        """
        # grab lengths of the sequences involved
        refLengths, e, em = backticks( 'fastalength %s' % self.sequence2 )
        refLength = int( refLengths[0].split()[0] )
        seqLengths, e, em = backticks( 'fastalength %s' % self.sequence1 )
        hitLengths = {}
        for line in seqLengths:
            values = line.split()
            hitLengths[ values[1] ] = int(values[0])
        #
        # TODO: handle multiple hits
        #
        print '%s %s' % ( os.path.realpath(self.sequence2), \
                          os.path.realpath(self.sequence1) )
        print 'NUCMER'
        for query_id, hitList in hits.iteritems():
            if len(hitList)==0: continue
            for hit in hitList:
                print hit.toDelta( refLength, hitLengths[hit.query_id] )

    def __setupMultipleHitPolicy(self):
        multiplePolicy = self.options.multiple.split(",")
        multiplePolicy[0] =  multiplePolicy[0].strip()
        if multiplePolicy[0] == "random":
            self.multipleHitPolicy = RandomMultipleHit()
        elif  multiplePolicy[0] == "all":
            self.multipleHitPolicy = AllMultipleHit()
        elif  multiplePolicy[0] == 'leftpositive' :
            if len(multiplePolicy)<2:
                deltaz = 3.0
            else:
                deltaz = float( multiplePolicy[1] )
            self.multipleHitPolicy = LeftPositiveMultipleHit( deltaz )
        elif multiplePolicy[0] == 'leftpositivelong':
            if len(multiplePolicy)<2:
                zth = 3.0
            else:
                zth = float( multiplePolicy[1] )
            self.multipleHitPolicy = LeftPositiveLongestHit( zth )
        elif multiplePolicy[0] == 'deltaz':
            if len(multiplePolicy)<2:
                deltaz = 2.0
            else:
                deltaz = float( multiplePolicy[1] )
            self.multipleHitPolicy = DeltaZHit( deltaz )
        elif multiplePolicy[0] =='bestscore':
            negate = False
            # depending on the selected algorithm "best scores" might be
            # negative or positive
            if self.options.algorithm == "blasr":
                negate = True
                self.multipleHitPolicy = BestScoreHit( negate=negate )
            else:
                # For an algorithm other than blasr,
                # if placeRepeatsRandomly is not set, always return the very
                # first hit reached with the best score; else return a random
                # hit with the best score
                if self.options.placeRepeatsRandomly == False :
                    self.multipleHitPolicy = BestScoreHit( negate=negate )
                else :
                    self.multipleHitPolicy = BestScoreHitPlaceRandomly( negate=negate )

        else:
            print >>sys.stderr, "Don't know multiple hit policy: %s" % self.options.multiple
            raise SystemExit

    def __loadAdapterAnnotations( self ):
        self._adapters = None
        if self.options.filterAdapterOnly:
            self._adapters = AdapterRanges( self.reference )

    def _getReferenceLength(self):
#        if not self.reference.is_pseudo_entry:

        #Why just return the first? What about multiseq fasta?

        return self.reference.contigs[0].length
#        sRefLength, errCode, errorMessage = backticks( 'fastalength %s' % self.sequence2 )
#
#        if errCode:
#            print >>sys.stderr, 'Error while running fastalength %s' % self.sequence2
#            print >>sys.stderr, errorMessage
#            raise SystemExit
#
#        return int(sRefLength[0].split()[0])

    def _getNumberReads(self):
        sNumberReads, errCode, errorMessage = backticks( "grep -c '>' %s" % self.sequence1 )
        if errCode:
            print >>sys.stderr, "Error while running grep -c '>' %s" % self.sequence1
            print >>sys.stderr, errorMessage
            sys.exit(1)
        return int(sNumberReads[0])

    def _setupLog(self, level=logging.DEBUG):
        """
        Setup global logger to dump to stdout.

        This is purposely set to be a bit verbose because the current code
        uses the global logger, not module level logging.
        """
        str_formatter = '[%(levelname)s] %(asctime)-15s [%(module)s %(funcName)s %(lineno)d] %(message)s'
        handler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter(str_formatter)
        handler.setFormatter(formatter)
        log.addHandler(handler)
        log.setLevel(level)


    def run(self):
        """
        Main Driver
        """
        started_at = time.time()

#        if self.options.debug:
#            level = logging.DEBUG
#        elif self.options.info:
#            level = logging.INFO
#        else:
#            level = logging.INFO
#
#        self._setupLog(level)
        log.info("Running {v} of {f} with algo -> {a}".format(v=CompareSequences.VERSION, f=self.__class__.__name__, a=self.options.algorithm))
#        log.info("Logging Level set to {d}".format(d=level))

        log.info("Processing Options")

        self.__setupMultipleHitPolicy()
        self.__loadAdapterAnnotations()

        # do quick checks of the reference length and
        # number of reads to tune the data parallelization
        refLength = self._getReferenceLength()
        #        numReads = self._getNumberReads()
        numReads = 0

        service = None
        if self.options.algorithm == "exonerate":
            service = ExonerateService()
        elif self.options.algorithm == "yass":
            service = YassService()
        elif self.options.algorithm == "blasr":
            service = BlasrService()
            service.setDebug( self.options.debug )

            # Make life easier on blasr and always use a temporary file.
            # Each thread in blasr prints to a separate file to avoid
            # jumbling the output, so create the base temporary file name here.
            service.useTemporaryFile = 1
            if self.options.debug == False:
                service.keepTemp = False
                service.temporaryFileName = self.algorithmParameters['-out'] = tempfile.mkstemp(dir=self.options.tmpDir)[1]
            else:
                service.keepTemp = True
                if self.options.manualTempFileName == "":
                    service.temporaryFileName = self.algorithmParameters['-out'] = tempfile.mkstemp(dir=self.options.tmpDir)[1]
                else:
                    service.temporaryFileName = self.algorithmParameters['-out'] = self.options.manualTempFileName
            if (self.options.algorithmOpts):
                for opt in self.options.algorithmOpts:
                    if (opt[0] == '-titleTable'):
                        self.algorithmParameters['-titleTable'] = opt[1]
            if ('-titleTable' not in self.algorithmParameters):
                self.algorithmParameters['-titleTable'] = tempfile.mkstemp(dir=self.options.tmpDir)[1]


            if (self.options.dryrun):
                service.dryRun = True
            else:
                service.dryRun = False

            service.nproc = 1
            if self.options.nproc > 1:
                self.algorithmParameters['-nproc'] = str(self.options.nproc)
                service.nproc = self.options.nproc

        elif self.options.algorithm == "lastz":
            service = LastzService()
        elif self.options.algorithm == "bwasw":
            service = BwaswService()
            self.algorithmParameters['-regionTable'] = self.options.regionTable

        else:
            msg = "Algorithms other than exonerate,yass,blasr,lastz and bwa are not yet supported"
            log.error(msg)
            raise SystemExit(msg)


        service.setQuery( self.sequence1 )

        service.setTarget( self.reference.sequenceFiles[0] )

        service.setParameters( self.algorithmParameters )

        if self.options.algorithm != "blasr":
            service.setNWorkers( self.options.nproc )

        if self.options.algorithm == "blasr":
            service.setCircularReference( self.options.circularReference )

        if self.options.showAlignment:
            service.setKeepAlignmentStrings( True )

        if refLength>LARGE_REFERENCE_LENGTH:
            # this logic should be revisited for large references
            # (the following only scales with the number of reads
            #  for now since exonerate seems to have a memory leak
            #  with respect to this quantity)
            # probably has an average readlength dependence too
            # (this logic seems reasonable for 8-core 16GB nodes
            #  when the average readlength is ~400bp)
            # Plan is to use true distributed processing
            # when we start to care more about huge-scale alignments
            multiplier = int(8.0 * float(numReads) / 500000.0)
            if multiplier < 1:
                multiplier = 1

            if (self.options.algorithm == "exonerate"):
                service.setWorkerMultiplier( multiplier )
        # TODO at some point it might be nice to remove as much exonerate specific stuff as possible

        if self.options.debug:
            self._logMemory( 'Before hits' )

        #
        # all the real work happens here
        #

        log.info("Running Alignment Service")
        service.run( self.options.tmpDir )
        log.info("Alignment Service Complete")

        hits = service.getHits()

        if self.options.debug:
            self._logMemory( 'After hits' )

        if self.options.delta:
            self._printDeltaOutput( hits )
        else:
            if self.options.h5fn and self.options.h5fn != "":
                log.info("Writing cmp.h5")
                self._saveAlignmentInHDF5(hits, service.getCommandLine(), mode=self.options.h5mode)
                log.info("Writing cmp.h5 Complete")

        if self.options.debug:
            self._logMemory( 'Exiting' )

        self._cleanUp(service)


        run_time = time.time() - started_at
        log.info("Completed in {s:.2f} sec ({m:.2f} min) {f}".format(s=run_time, m=run_time/60.0, f=__file__))


    def _cleanUp(self,service):
        '''
        Delete tmp files
        '''
        if self.options.pls2fasta == True:
            os.remove( self.tmpFastaFileName )
        if not self.options.debug:
            if self.options.algorithm == "blasr":
                os.remove(service.temporaryFileName)
            if self._tmpRepos is not None:
                backticks( 'rm -rf %s' % self._tmpRepos )



def pct( num, den ):
    if den==0:
        return 0.0
    return float(num)/float(den) * 100.0

def calculateZ( zcalculator, hit ):
    if not zcalculator:
        return 0.0
    if (not hit.hasMatchInfo and hit.hasAlignments):
        raise SystemExit, "SW error, hit needs alignment and match info to calc Z-score"
    al = abs( hit.query_end - hit.query_start )
    nCorrect = al - hit.nIns - hit.nDel - hit.nMismatch
    if (nCorrect < 0):
        nCorrect = 0

    try:
        if zcalculator.use_ref_correction:
            #seq = hit.getQueryAlignment( read ).replace('-','')
            seq = hit.alignedQuery.replace('-','')
            h = calculateEntropyForZscore( seq )
            ref_min = min(hit.target_start,hit.target_end)
            ref_max = max(hit.target_start,hit.target_end)
            z = zcalculator.calcZ( al, pct(nCorrect,al), entropy=h, \
                ref_min=ref_min, ref_max=ref_max )
        elif zcalculator.use_h_correction:
            #seq = hit.getQueryAlignment( read ).replace('-','')
            seq = hit.alignedQuery.replace('-','')
            h = calculateEntropyForZscore( seq )
            #print >>sys.stderr, 'seq = %s' % seq
            #print >>sys.stderr, 'h = %f' % h
            z = zcalculator.calcZ( al, pct(nCorrect,al), entropy=h )
        else:
            z = zcalculator.calcZ( al, pct(nCorrect,al) )

    except Exception, e:
        # shouldn't happen but there can be system instabilities (!!)
        z = -10.0

    return z

class IMultipleHitPolicy:
    """Interface for business logic which filters a list of candidate hits."""
    def handleMultipleHits(self, hitList ):
        """Implementors supply logic for pruning a list of multiple candidate hits.
        Returns a list of pruned hits.
        """
        pass

    def setZCalculator(self, zcalculator):
        pass

class RandomMultipleHit(IMultipleHitPolicy):
    """Multiple hit policy which returns a randomly selected single hit."""
    def handleMultipleHits(self, hitList ):
        if len(hitList)==1: return hitList
        i = random.randint( 0, len(hitList)-1 )
        return [ hitList[i] ]

class LeftPositiveMultipleHit(IMultipleHitPolicy):
    """Multiple hit policy which selects the left most positive strand hit
    as long as it is within a prescribed Z-threshold of the top hit."""
    def __init__(self, deltaz):
        self.deltaz = deltaz
        self.zcalculator = None

    def setZCalculator(self, zcalculator):
        self.zcalculator = zcalculator

    def handleMultipleHits(self, hitList ):
        """Calculates Z-score for each hit and only returns a hit
        if all alternative hits are more than self.deltaz away."""
        #print >>sys.stderr, "considering %d hits" % len(hitList)
        zs = [ self._calcZ(hit) for hit in hitList ]
        if len(zs):
            zmax = max(zs)
        else:
            # We removed zscore from the file.
            zmax = 30.0
        zmax = max(zs)
        maxcount = 0
        deltamin = 1.0e10
        bestHit = None
        leftHit = None
        validHits = 0
        for z, hit in zip( zs, hitList ):
            if ((leftHit == None or hit.target_start < leftHit.target_start) and \
                    hit.target_strand == '+' and hit.query_strand == '+' and \
                    zmax - z < self.deltaz):
                    bestHit = hit
                    leftHit = hit
                    validHits += 1
        if (validHits > 0):
            return [ bestHit ]
        return []

    def _calcZ(self, hit):
        return calculateZ( self.zcalculator, hit )

class LeftPositiveLongestHit(IMultipleHitPolicy):
    """Multiple hit policy which selects the longest postive strand alignemnt with Z>3.
    Where there is degeneracy, the left most one is picked """
    def __init__(self, zth=None):
        if zth != None:
            self.zth = zth
        else:
            self.zth = 3
        self.zcalculator = None

    def setZCalculator(self, zcalculator):
        self.zcalculator = zcalculator

    def handleMultipleHits(self, hitList ):
        """Calculates Z-score for each hit and only returns a hit
        if all alternative hits are more than self.deltaz away."""
        #print >>sys.stderr, "considering %d hits" % len(hitList)
        zs = [ self._calcZ(hit) for hit in hitList ]
        tlen = [ hit.target_end - hit.target_start for hit in hitList ]
        ts = [ hit.target_start for hit in hitList ]
        zmax = max(zs)
        filteredHit = zip( tlen, ts, zs, hitList )
        filteredHit = [ x for x in filteredHit if x[2] > self.zth and x[3].target_strand == '+' and x[3].query_strand == '+' ]
        if len(filteredHit) != 0:
            filteredHit.sort(key=lambda x:(-x[0],x[1]) ) #sort by length then by position
            return [filteredHit[0][3]]
        else:
            return []

    def _calcZ(self, hit):
        return calculateZ( self.zcalculator, hit )

class AllMultipleHit(IMultipleHitPolicy):
    """Multiple hit policy which allows all hits to be reported."""
    def handleMultipleHits(self, hitList ):
        return hitList

class DeltaZHit( IMultipleHitPolicy ):
    """Multiple hit policy which requires alternate hits to be
    worse by more than a prescribed Z-threshold in order to report a hit."""
    def __init__(self, deltaz):
        self.deltaz = deltaz
        self.zcalculator = None

    def setZCalculator(self, zcalculator):
        self.zcalculator = zcalculator

    def handleMultipleHits(self, hitList ):
        """Calculates Z-score for each hit and only returns a hit
        if all alternative hits are more than self.deltaz away."""
        #print >>sys.stderr, "considering %d hits" % len(hitList)
        zs = [ self._calcZ(hit) for hit in hitList ]
        zmax = max(zs)
        maxcount = 0
        deltamin = 1.0e10
        bestHit = None
        for z, hit in zip( zs, hitList ):
            if z==zmax:
                maxcount += 1
                bestHit = hit
                continue
            if zmax-z < deltamin:
                deltamin = zmax-z
        if maxcount>1:
            return []
        if deltamin < self.deltaz:
            return []
        return [ bestHit ]

    def _calcZ(self, hit):
        return calculateZ( self.zcalculator, hit )

class BestScoreHit( IMultipleHitPolicy ):
    """Multiple hit policy which reports the very first hit
       with the best score."""
    def __init__(self, negate=False):
        self._sign = -1.0 if negate else 1.0

    def handleMultipleHits( self, hitList ):
        if len(hitList)<=1: return hitList
        maxScore = -1.0e10
        iMax = 0
        for i,h in enumerate(hitList):
            s = h.score * self._sign
            if s > maxScore:
                iMax = i
                maxScore = s
        return [ hitList[iMax] ]

class BestScoreHitPlaceRandomly( IMultipleHitPolicy ):
    """Multiple hit policy which reports a random hit with the
       the best score."""
    def __init__(self, negate=False):
        self._sign = -1.0 if negate else 1.0

    def handleMultipleHits( self, hitList ):
        if len(hitList) <= 1: return hitList
        maxScore = -1.0e10
        hitIndexs = []
        hitList.sort(key=lambda x:(x.query_start, x.query_end, x.target_start, x.target_end))
        for i,h in enumerate(hitList):
            s = h.score * self._sign
            if s > maxScore:
                maxScore = s
                hitIndexs = [i]
            elif s ==  maxScore:
                hitIndexs.append(i)
        iRand = random.randint( 0, len(hitIndexs)-1 )
        return [ hitList[ hitIndexs[iRand] ] ]


class AdapterRanges:
    def __init__( self, reference ):
        self._ranges = Ranges()
        if len(reference.annotations)!=1:
            return
        gffFile = reference.annotations[0]._file
        self._loadFromGff( gffFile )

    def _loadFromGff( self, gffFile ):
        reader = GffReader( gffFile )
        for record in reader:
            if record.type=='adapter':
                # convert to 0-based, exclusive
                r = Range( record.start-1, record.end )
                self._ranges.addRange(r)
        reader.close()

    def checkAdapterOnly( self, hit ):
        """Returns True if this hit overlaps an adapter region and
        if the combined region is no bigger than
        adapter region + FUZZY_OVERLAP"""
        FUZZY_OVERLAP = 20
        hitRange = Range( hit.target_start, hit.target_end )
        for r in self._ranges:
            if r.intersects( hitRange ):
                lengthUnion = max(r.end,hitRange.end) - min(r.start,hitRange.start)
                if lengthUnion < len(r)+FUZZY_OVERLAP:
                    return True
        return False

class HomopolymerGroup:
    """Used by the randomize deletions algorithm."""
    def __init__( self, qCh, tCh ):
        self.ref = tCh
        self.n = 1
        self.bases = [ qCh ]

    def add( self, qCh ):
        self.bases.append( qCh )
        self.n += 1

    def shuffle( self, rObj ):
        """randomizes the location of gaps, but keeps the
        underlying query sequence the same"""
        if self.n==1:
            return self.bases[0], self.ref
        if '-' not in self.bases:
            return ''.join(self.bases), self.n * self.ref
        calls = ''.join([ch for ch in self.bases if ch!='-'])
        nGaps = self.n - len(calls)
        newBases = [ 'x' for i in xrange(self.n) ]
        for i in rObj.sample( xrange(self.n), nGaps ):
            newBases[i] = '-'
        j = 0
        for i in xrange(self.n):
            if newBases[i]=='x':
                newBases[i] = calls[j]
                j += 1
        return ''.join(newBases), self.n*self.ref

    def __str__( self ):
        """for debugging"""
        return '%s%d:%s' % ( self.ref, self.n, "".join(self.bases) )

def randomize_deletions( q, t ):
    """For the SMRT sequencing error model it can make sense
    to randomize the location of deletions within a homopolymer
    stretch.  For aligned query q and target t, return a version
    where the location of gaps in homopolymer stretches of t
    is randomized.  The original order of basecalls in q is
    preserved."""
    r = random.Random()
    group = None
    groups = []
    for ch1, ch2 in zip(q,t):
        if not group:
            group = HomopolymerGroup( ch1, ch2 )
            continue
        if ch2!=group.ref:
            groups.append(group)
            group = HomopolymerGroup( ch1, ch2 )
        else:
            group.add( ch1 )
    if group:
        groups.append(group)

    a = []
    for g in groups:
        a.append( g.shuffle(r) )
    return [ ''.join(s) for s in zip(*a) ]

if __name__=='__main__':
    app = CompareSequences( sys.argv )
    if app.options.profile:
        import cProfile
        cProfile.run( 'app.run()', 'compareSequences.profile' )
    else:
        app.run()
