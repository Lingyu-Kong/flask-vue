<template>
  <div>
    <p>Home page</p>
    <!-- <p>Random number from backend: {{ randomNumber }}</p>
    <button @click="getRandom">New random number</button> -->
    <!-- 上传文件 -->
    <input type="file" @change="onFileChange" />
    <button @click="uploadFile">上传文件</button>
    <!-- 展示计算结果 -->
    <div v-if="result !== null">
        <h3>计算结果：</h3>
        <pre>{{ result }}</pre>
    </div>
  </div>
</template>

<script>
import axios from 'axios'
export default {
  data () {
    return {
      randomNumber: 0,
      selectedFile: null,
      result: null
    }
  },
  methods: {
    getRandom () {
      this.randomNumber = this.getRandomFromBackend()
    },
    getRandomFromBackend () {
      const path = `http://localhost:5000/api/random`
      axios
        .get(path)
        .then((response) => {
          this.randomNumber = response.data.randomNumber
        })
        .catch((error) => {
          console.log(error)
        })
    },
    onFileChange (event) {
      this.selectedFile = event.target.files[0]
    },
    async uploadFile () {
      if (!this.selectedFile) {
        alert('请先选择一个文件')
        return
      }
      const formData = new FormData()
      formData.append('file', this.selectedFile)
      try {
        const response = await axios.post('http://localhost:5000/api/upload', formData, {
          headers: {
            'Content-Type': 'multipart/form-data'
          }
        })
        this.result = response.data
      } catch (error) {
        console.error('Error uploading file:', error)
      }
    }
  },
  created () {
    this.getRandom()
  }
}
</script>
