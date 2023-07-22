<template>
  <div>
    <h2>Home Page</h2>
    <input type='file' @change='onFileChange' />
    <button @click='upload'>Upload File</button>
    <button v-if="downloadLink" @click="() => downloadResults(downloadLink, result_file)">Download Results</button>
    <div v-if='waiting'>Processing...</div>
    <!-- 展示计算结果 -->
    <div v-if="result !== null" class="centered-container">
      <!-- <h1>Single Row Table</h1>   -->
      <ResultsTable :dictionary="result" />
    </div>
    <div v-if="result !== null">
      <label for="steps">Steps:</label>
      <input type="number" id="steps" v-model.number="steps" />
      <label for="fmax">Fmax:</label>
      <input type="number" id="fmax" step="0.01" v-model.number="fmax" />
      <button @click="upload_relax">Perform Relax</button>
      <button v-if="downloadLink_relax" @click="() => downloadResults(downloadLink_relax, relax_file)">Download Trajectory</button>
    </div>
    <div v-if='waiting_relax'>Relaxing...</div>
  </div>
</template>

<style>
.centered-container {
  display: flex;
  justify-content: center;
  align-items: center;
  min-height: 100%; /* Optional: Set a minimum height if needed */
}
</style>

<script>

import axios from 'axios'
import ResultsTable from '@/components/ResultsTable.vue'

export default {
  components: {
    ResultsTable
  },
  data () {
    return {
      file: null,
      waiting: false,
      downloadLink: '',
      result: null,
      result_file: '',
      waiting_relax: false,
      downloadLink_relax: '',
      relax_file: '',
      steps: 100,
      fmax: 0.1
    }
  },
  methods: {
    onFileChange (e) {
      this.file = e.target.files[0]
      this.downloadLink = ''
      this.result = null
      this.result_file = ''
      this.downloadLink_relax = ''
      this.relax_file = ''
    },
    async upload () {
      if (!this.file) {
        alert('Please select a file for upload')
        return
      }

      const formData = new FormData()
      formData.append('file', this.file)

      try {
        this.waiting = true
        const response = await fetch('http://localhost:5000/upload', {
          method: 'POST',
          body: formData
        })
        const data = await response.json()
        if (data.status === 'success') {
          this.downloadLink = `http://localhost:5000/download/${data.filename}`
          this.result = data.result
          this.result_file = data.filename
        } else {
          alert(data.message)
        }
      } catch (error) {
        alert('Error in uploading file')
      } finally {
        this.waiting = false
      }
    },
    async downloadResults (downloadLink, fileName) {
      try {
        const response = await axios.get(downloadLink, { responseType: 'blob' })
        const url = window.URL.createObjectURL(new Blob([response.data]))
        const link = document.createElement('a')
        link.href = url
        link.setAttribute('download', fileName) // Use the original file name for the downloaded file
        document.body.appendChild(link)
        link.click()
        link.remove()
      } catch (error) {
        alert('Error downloading the file')
      }
    },
    async upload_relax () {
      if (!this.file) {
        alert('Please select a file for upload')
        return
      }

      const formData = new FormData()
      formData.append('file', this.file)
      formData.append('steps', this.steps)
      formData.append('fmax', this.fmax)

      try {
        this.waiting_relax = true
        const response = await fetch('http://localhost:5000/upload_relax', {
          method: 'POST',
          body: formData
        })
        const data = await response.json()
        if (data.status === 'success') {
          this.downloadLink_relax = `http://localhost:5000/download/${data.filename}`
          this.relax_file = data.filename
        } else {
          alert(data.message)
        }
      } catch (error) {
        alert('Error in uploading file')
      } finally {
        this.waiting_relax = false
      }
    }
  }
}
</script>
