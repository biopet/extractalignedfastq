package nl.biopet.tools.extractalignedfastq

import java.io.File

case class Args(inputBam: File = new File(""),
                intervals: List[String] = List.empty[String],
                inputFastq1: File = new File(""),
                inputFastq2: Option[File] = None,
                outputFastq1: File = new File(""),
                outputFastq2: Option[File] = None,
                minMapQ: Int = 0,
                commonSuffixLength: Int = 0)
