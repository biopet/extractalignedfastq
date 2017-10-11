package nl.biopet.tools.extractalignedfastq

import nl.biopet.test.BiopetTest
import org.testng.annotations.Test

class ExtractAlignedFastqTest extends BiopetTest {
  @Test
  def testNoArgs(): Unit = {
    intercept[IllegalArgumentException] {
      ExtractAlignedFastq.main(Array())
    }
  }
}
