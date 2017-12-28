package nl.biopet.tools.extractalignedfastq

import java.io.File

import htsjdk.samtools.{QueryInterval, SamReaderFactory, ValidationStringency}
import htsjdk.samtools.fastq.{BasicFastqWriter, FastqReader, FastqRecord}
import htsjdk.samtools.util.Interval
import nl.biopet.utils.tool.ToolCommand

import scala.collection.mutable.{Set => MSet}
import scala.collection.JavaConverters._

object ExtractAlignedFastq extends ToolCommand[Args] {
  def emptyArgs: Args = Args()
  def argsParser = new ArgsParser(this)
  def main(args: Array[String]): Unit = {
    val cmdArgs = cmdArrayToArgs(args)

    // Require input files to exist
    require(cmdArgs.inputBam.exists(),"Input BAM file not found")
    require(cmdArgs.inputFastq1.exists(),"Input FASTQ file 1 not found")
    if (cmdArgs.inputFastq2.isDefined){
      require(cmdArgs.inputFastq2.get.exists(),"Input FASTQ file 2 not found")
    }

    logger.info("Start")

    val memFunc = makeMembershipFunction(
      iv = makeIntervalFromString(cmdArgs.intervals),
      inAln = cmdArgs.inputBam,
      minMapQ = cmdArgs.minMapQ,
      commonSuffixLength = cmdArgs.commonSuffixLength
    )

    logger.info("Writing to output file(s) ...")
    (cmdArgs.inputFastq2, cmdArgs.outputFastq2) match {

      case (None, None) =>
        val in = new FastqReader(cmdArgs.inputFastq1)
        val out = new BasicFastqWriter(cmdArgs.outputFastq1)
        extractReads(memFunc, in, out)
        in.close()
        out.close()

      case (Some(i2), Some(o2)) =>
        val in1 = new FastqReader(cmdArgs.inputFastq1)
        val in2 = new FastqReader(i2)
        val out1 = new BasicFastqWriter(cmdArgs.outputFastq1)
        val out2 = new BasicFastqWriter(o2)
        extractReads(memFunc, in1, out1, in2, out2)
        in1.close()
        in2.close()
        out1.close()
        out2.close()

      case _ => ; // handled by the command line config check above
    }

    logger.info("Done")
  }

  /** type alias for Fastq input (may or may not be paired) */
  type FastqInput = (FastqRecord, Option[FastqRecord])

  /** Get the FastqRecord ID */
  def fastqId(rec: FastqRecord): String = rec.getReadName.split(" ")(0)

  /**
    * Function to create iterator over Interval given input interval string
    *
    * Valid interval strings are either of these:
    *    - chr5:10000-11000
    *    - chr5:10,000-11,000
    *    - chr5:10.000-11.000
    *    - chr5:10000-11,000
    * In all cases above, the region span base #10,000 to base number #11,000 in chromosome 5
    * (first base is numbered 1)
    *
    * An interval string with a single base is also allowed:
    *    - chr5:10000
    *    - chr5:10,000
    *    - chr5:10.000
    *
    * @param inStrings iterable yielding input interval string
    */
  def makeIntervalFromString(inStrings: Iterable[String]): Iterator[Interval] = {

    // FIXME: can we combine these two patterns into one regex?
    // matches intervals with start and end coordinates
    val ptn1 = """([\w_-]+):([\d.,]+)-([\d.,]+)""".r
    // matches intervals with start coordinate only
    val ptn2 = """([\w_-]+):([\d.,]+)""".r
    // make ints from coordinate strings
    // NOTE: while it is possible for coordinates to exceed Int.MaxValue, we are limited
    // by the Interval constructor only accepting ints
    def intFromCoord(s: String): Int =
      s.replaceAll(",", "").replaceAll("\\.", "").toInt

    inStrings.map {
      case ptn1(chr, start, end) if intFromCoord(end) >= intFromCoord(start) =>
        new Interval(chr, intFromCoord(start), intFromCoord(end))
      case ptn2(chr, start) =>
        val startCoord = intFromCoord(start)
        new Interval(chr, startCoord, startCoord)
      case otherwise =>
        throw new IllegalArgumentException(
          "Invalid interval string: " + otherwise)
    }.toIterator
  }

  /**
    * Function to create object that checks whether a given FASTQ record is mapped
    * to the given interval or not
    *
    * @param iv iterable yielding features to check
    * @param inAln input SAM/BAM file
    * @param minMapQ minimum mapping quality of read to include
    * @param commonSuffixLength length of suffix common to all read pairs
    * @return
    */
  def makeMembershipFunction(
      iv: Iterator[Interval],
      inAln: File,
      minMapQ: Int = 0,
      commonSuffixLength: Int = 0): (FastqInput => Boolean) = {

    val inAlnReader = SamReaderFactory
      .make()
      .validationStringency(ValidationStringency.LENIENT)
      .open(inAln)
    require(inAlnReader.hasIndex)

    def getSequenceIndex(name: String): Int =
      inAlnReader.getFileHeader.getSequenceIndex(name) match {
        case x if x >= 0 =>
          x
        case _ =>
          throw new IllegalArgumentException(
            "Chromosome " + name + " is not found in the alignment file")
      }

    val queries: Array[QueryInterval] = iv.toList
    // transform to QueryInterval
      .map(x =>
        new QueryInterval(getSequenceIndex(x.getContig), x.getStart, x.getEnd))
      // sort Interval
      .sortBy(x => (x.referenceIndex, x.start, x.end))
      // cast to array
      .toArray

    lazy val selected: MSet[String] = inAlnReader
    // query BAM file for overlapping reads
      .queryOverlapping(queries)
      // for Scala compatibility
      .asScala
      // filter based on mapping quality
      .filter(x => x.getMappingQuality >= minMapQ)
      // iteratively add read name to the selected set
      .foldLeft(MSet.empty[String])(
        (acc, x) => {
          logger.debug("Adding " + x.getReadName + " to set ...")
          acc += x.getReadName
        }
      )

    (pair: FastqInput) =>
      pair._2 match {
        case None => selected.contains(fastqId(pair._1))
        case Some(x) =>
          val rec1Id = fastqId(pair._1)
          require(commonSuffixLength < rec1Id.length)
          require(commonSuffixLength < fastqId(x).length)
          selected.contains(rec1Id.dropRight(commonSuffixLength))
      }
  }

  /**
    * Extracts reads from the given input Fastq file and writes to a new output Fastq file
    *
    * @param memFunc Predicate for extracting reads. If evaluates to true, the read is extracted.
    * @param inputFastq1 Input [[FastqReader]] object.
    * @param outputFastq1 Output [[BasicFastqWriter]] object.
    */
  def extractReads(memFunc: FastqInput => Boolean,
                   inputFastq1: FastqReader,
                   outputFastq1: BasicFastqWriter): Unit =
    inputFastq1.iterator.asScala
      .filter(rec => memFunc((rec, None)))
      .foreach(rec => outputFastq1.write(rec))

  /**
    * Extracts reads from the given input Fastq pairs and writes to new output Fastq pair files
    *
    * @param memFunc Predicate for extracting reads. If evaluates to true, the read is extracted.
    * @param inputFastq1 Input [[FastqReader]] object for pair 1.
    * @param outputFastq1 Input [[FastqReader]] object for pair 2.
    * @param inputFastq2 Output [[BasicFastqWriter]] object for pair 1.
    * @param outputFastq2 Output [[BasicFastqWriter]] object for pair 2.
    */
  def extractReads(memFunc: FastqInput => Boolean,
                   inputFastq1: FastqReader,
                   outputFastq1: BasicFastqWriter,
                   inputFastq2: FastqReader,
                   outputFastq2: BasicFastqWriter): Unit =
    inputFastq1.iterator.asScala
      .zip(inputFastq2.iterator.asScala)
      .filter(rec => memFunc(rec._1, Some(rec._2)))
      .foreach(rec => {
        outputFastq1.write(rec._1)
        outputFastq2.write(rec._2)
      })

  def descriptionText: String =
    s"""
       |This tool extracts reads from a BAM file based on alignment intervals.
       |E.g if one is interested in a specific location
       |this tool extracts the full reads from the location.
       |The tool is also very useful to create test data sets.
     """.stripMargin

  def manualText: String =
    """
      |This tool creates FASTQ file(s) containing reads mapped to the given alignment intervals. A set of FASTQ files that was
      |used in creating the BAM file is also required since this is used for retrieving full sequences of FASTQ records which
      |map to the given region. This is useful since some of the records may have undergone modifications such as quality
      |trimming before alignment. In this case, retrieving the aligned SAM records will only give the modified sequence.
      |
      |* Note that this tool works for single end and paired end data. The below example can be easily extended for
      |paired end data. The only thing one should add is: --in2 myFastq_R2.fastq --out2 myOutFastq_R2.fastq
      |* The interval is just a genomic position or multiple genomic positions wherefrom one wants to extract the reads.
      |
    """.stripMargin

  def exampleText: String =
    s"""
       |To extract reads from `myBam.bam` that originate from `myFastq_r1.fasta` and align to chr5 at positions 100-200:
       |${example("--input_file",
                  "myBam.bam",
                  "--in1",
                  "myFastq_R1.fastq",
                  "--interval",
                  "chr5:100-200")}
     """.stripMargin
}
